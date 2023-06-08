import sys
import os
import pandas as pd
import numpy as np
from tqdm import tqdm

import h5py
import pyfaidx
from gtfparse import read_gtf


def co_to_idx(start, end, strand):
    return start - 1, end


def slice_gen(
    seq,
    start,
    end,
    strand,
    co=True,
    to_vec=True,
    seq_dict={"A": 0, "T": 1, "C": 2, "G": 3, "N": 4},
    comp_dict={0: 1, 1: 0, 2: 3, 3: 2, 4: 4},
):
    if co:
        start, end = co_to_idx(start, end, strand)
    sl = seq[start:end].seq

    if to_vec:
        sl = list(map(lambda x: seq_dict[x], sl))

    if strand in ["-", -1, False]:
        if comp_dict is not None:
            sl = list(map(lambda x: comp_dict[x], sl))[::-1]
        else:
            sl = sl[::-1]

    return np.array(sl)


def parse_transcriptome(gtf, fasta, chromsizes, h5_out, meta_dir):
    genome = pyfaidx.Fasta(fasta)
    contig_list = pd.read_csv(chromsizes, sep="\t", header=None, dtype={0: str})[0]
    gtf = read_gtf(gtf)

    print("Extracting transcripts and metadata...")
    os.makedirs(os.path.join(meta_dir, "numpy"), exist_ok=True)
    os.makedirs(os.path.join(meta_dir, "metadata"), exist_ok=True)

    for contig in contig_list:
        print(f"{contig}...")
        tran_ids = []
        samples = []
        poi_contig = []
        gtf_set = gtf[gtf.seqname == str(contig)]
        tr_set = gtf_set["transcript_id"].unique()[1:]

        tr_datas = []
        headers = [
            "transcript_id",
            "gene_id",
            "gene_name",
            "strand",
            "transcript_biotype",
            "tag",
            "transcript_support_level",
        ]
        for i, tr_idx in tqdm(enumerate(tr_set), total=len(tr_set)):
            # obtain transcript information
            gtf_tr = gtf_set[gtf_set.transcript_id == tr_idx]
            tr_data = gtf_tr.loc[gtf_tr.feature == "transcript", headers].iloc[0]

            # obtain and sort exon information
            exons = gtf_tr[gtf_tr.feature == "exon"].copy()
            exons.exon_number = exons.exon_number.astype(int)
            exons = exons.sort_values("exon_number")

            # obtain TISs
            TISs = gtf_tr[gtf_tr.feature == "start_codon"].sort_values("exon_number")
            TISs = TISs.iloc[0:1]

            # obtain TTSs
            CDSs = gtf_tr[gtf_tr.feature == "CDS"].sort_values("exon_number")

            exon_seqs = []
            target_seqs = []
            poi_tr = []
            tr_len = 0

            for j, exon in exons.iterrows():
                exon_seq = slice_gen(
                    genome[contig], exon.start, exon.end, exon.strand, to_vec=True
                ).astype(np.int16)
                if exon.strand == "+":
                    exon_coords = [exon.start, exon.end]
                else:
                    exon_coords = [exon.end, exon.start]
                poi_tr.append(["exon", tr_len, tr_len + len(exon_seq)] + exon_coords)
                exon_seqs.append(exon_seq)
                target_seq = np.full(exon_seq.shape, False)

                if len(TISs) > 0:
                    for k, tis in TISs[
                        np.logical_and(TISs.start >= exon.start, TISs.start <= exon.end)
                    ].iterrows():
                        # complexity related to varying strands
                        if tis.strand == "+":
                            target_seq[tis.start - exon.start] = 1
                            tis_idx = tr_len + tis.start - exon.start
                            tis_coords = [tis.start, tis.start + 1]

                        else:
                            target_seq[exon.end - tis.end] = 1
                            tis_idx = tr_len + exon.end - tis.end
                            tis_coords = [tis.end, tis.end - 1]
                        tr_data["tis_idx"] = tis_idx
                        idxs = np.hstack(([tis_idx, tis_idx + 1], tis_coords))
                        poi_tr.append(["tis"] + list(idxs))
                        mask = np.logical_and(
                            np.logical_or(
                                gtf_tr.start == tis.start, gtf_tr.end == tis.end
                            ),
                            gtf_tr.feature == "CDS",
                        )
                        if mask.sum() > 0:
                            tr_data["prot_ID"] = gtf_tr.loc[mask, "protein_id"].item()
                target_seqs.append(target_seq)
                tr_len += len(target_seq)

            if len(CDSs) > 0:
                cds_length = 0
                for k, cds in CDSs.iterrows():
                    # complexity related to varying strands
                    cds_length += abs(cds.start - cds.end) + 1
                if cds.strand == "+":
                    poi_tr.append(["tts"] + list(idxs + cds_length))
                else:
                    poi_tr.append(
                        ["tts"]
                        + list(idxs[:2] + cds_length)
                        + list(idxs[2:] - cds_length)
                    )
                tr_data["tts_idx"] = idxs[0] + cds_length

            tr_datas.append(tr_data)
            sample = np.vstack((np.concatenate(exon_seqs), np.concatenate(target_seqs)))

            samples.append(sample)
            tran_ids.append(tr_idx)
            poi_contig.append(poi_tr)

        final = np.array([tran_ids, samples, poi_contig], dtype=object).T
        if len(tr_datas) == 0:
            metadata = pd.DataFrame(columns=headers)
        else:
            metadata = pd.concat(tr_datas, axis=1).T
        metadata.to_csv(os.path.join(meta_dir, "metadata", f"{contig}.csv"))
        np.save(os.path.join(meta_dir, "numpy", f"{contig}.npy"), final)

    print("Save data in hdf5 files...")

    ## Create h5py database
    dt = h5py.vlen_dtype(np.dtype("int8"))
    contig_type = f"<S{contig_list.str.len().max()}"

    f = h5py.File(f"{h5_out}.h5", "w")
    if "transcript" not in f.keys():
        grp = f.create_group("transcript")
    else:
        grp = f["transcript"]
    if "metadata" not in f["transcript"]:
        f["transcript"].create_group("metadata")

    seqs = []
    tiss = []
    contigs = []
    ids = []
    len_trs = []
    for contig in contig_list:
        print(f"{contig}...")
        x = np.load(os.path.join(meta_dir, "numpy", f"{contig}.npy"), allow_pickle=True)
        seqs += [t[0] for t in x[:, 1]]
        tiss += [t[1] for t in x[:, 1]]
        contigs.append(np.full(len(x), contig, dtype=contig_type))
        ids.append(x[:, 0].astype("S15"))

    grp.create_dataset("seq", data=seqs, dtype=dt)
    grp.create_dataset("tis", data=tiss, dtype=dt)
    grp.create_dataset("contig", data=np.hstack(contigs))
    grp.create_dataset("id", data=np.hstack(ids))

    len_trs += [len(seq) for seq in seqs]
    grp.create_dataset("tr_len", data=len_trs)

    f.close()


if __name__ == "__main__":
    if len(sys.argv) < 6 or len(sys.argv) > 6:
        print("process_data [gtf] [fasta] [chrom.sizes], [h5_out] [meta_dir]")
    else:
        parse_transcriptome(*sys.argv[1:])
