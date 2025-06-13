#!/bin/bash

# Specify the directory containing the .sbatch scripts
script_directory=$1

# Get the list of sbatch scripts in the directory, sorted numerically
sbatch_scripts=($(ls -v "$script_directory"/*.sbatch))

cp archive_job_output.sh ${script_directory}/archive_job_output.sh

# Start the first 10 jobs without waiting
for ((i=0; i<10 && i<${#sbatch_scripts[@]}; i++)); do
    job_id=$(sbatch "${sbatch_scripts[i]}" | awk '{print $4}')
    output_directory=$(scontrol show job $job_id | awk -F= '/WorkDir/{print $2}')

    # Submit the archive job for the completed job
    tar_job_id=$(sbatch --dependency=afterok:$job_id "$script_directory/archive_job_output.sh" "$output_directory" | awk '{print $4}')
    tar_job_ids+=($tar_job_id)
done

# Start additional jobs as the previous ones finish
for ((i=10; i<${#sbatch_scripts[@]}; i++)); do
    # Wait for the first job to finish and be archive
    wait_tar_job_id=${tar_job_ids[0]}
    tar_job_ids=("${tar_job_ids[@]:1}")

    # Start the new job without waiting
    new_job_id=$(sbatch --dependency=afterok:$wait_tar_job_id "${sbatch_scripts[i]}" | awk '{print $4}')
    new_output_directory=$(scontrol show job $new_job_id | awk -F= '/WorkDir/{print $2}')
    
    # Submit the archive job for the completed job
    tar_job_id=$(sbatch --dependency=afterok:$new_job_id "$script_directory/archive_job_output.sh" "$new_output_directory" | awk '{print $4}')
    tar_job_ids+=($tar_job_id)
done

