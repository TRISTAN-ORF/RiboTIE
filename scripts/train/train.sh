#!/bin/bash

conda activate riboformer_env

transcript_transformer train 'template.json' --val 2 14 --test 1 7 13 19 --gpus 1 --max_epochs 40 --name 'template_fold_1' --dim 42 --depth 6 --heads 6 --dim_head 16 --local_attn_heads 4 --log_dir ../models/temp_fold_1
transcript_transformer train 'template.json' --val 1 13 --test 2 8 14 20 --gpus 1 --max_epochs 40 --name 'template_fold_2' --dim 42 --depth 6 --heads 6 --dim_head 16 --local_attn_heads 4 --log_dir ../models/temp_fold_2
transcript_transformer train 'template.json' --val 1 13 --test 3 9 15 21 --gpus 1 --max_epochs 40 --name 'template_fold_3' --dim 42 --depth 6 --heads 6 --dim_head 16 --local_attn_heads 4 --log_dir ../models/temp_fold_3
transcript_transformer train 'template.json' --val 1 13 --test 4 10 16 22 --gpus 1 --max_epochs 40 --name 'template_fold_4' --dim 42 --depth 6 --heads 6 --dim_head 16 --local_attn_heads 4 --log_dir ../models/temp_fold_4
transcript_transformer train 'template.json' --val 1 13 --test 5 11 17 X --gpus 1 --max_epochs 40 --name 'template_fold_5' --dim 42 --depth 6 --heads 6 --dim_head 16 --local_attn_heads 4 --log_dir ../models/temp_fold_5
transcript_transformer train 'template.json' --val 1 13 --test 6 12 18 Y --gpus 1 --max_epochs 40 --name 'template_fold_6' --dim 42 --depth 6 --heads 6 --dim_head 16 --local_attn_heads 4 --log_dir ../models/temp_fold_6