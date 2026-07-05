import os, sys, json
import torch
import numpy as np
import subprocess


dev_count = torch.cuda.device_count()
# del best_dev_id
print(f"CUDA device count: {dev_count}")
best_dev_id = 0
if dev_count > 0:
    free_fracs = []
    candidate_ids = []
    for dev_id in range(dev_count):

        try:
            m_free, m_total = torch.cuda.mem_get_info(dev_id)
            free_fracs.append(m_free / max(m_total, 1))
            candidate_ids.append(dev_id)
        except Exception:
            free_fracs.append(-1.0)
            candidate_ids.append(dev_id)
    if free_fracs:
        cand_idx = int(np.argmax(free_fracs))
        if free_fracs[cand_idx] >= 0:
            best_dev_id = candidate_ids[cand_idx]
# print gpu usage with PIDs and users
try:
    gpu_result = subprocess.run(
        [
            'nvidia-smi',
            '--query-gpu=index,memory.free,memory.total,utilization.gpu,utilization.memory',
            '--format=csv,noheader,nounits'
        ],
        stdout=subprocess.PIPE,
        text=True,
        check=False,
    )
    print("GPU Usage:")
    print(gpu_result.stdout)

    proc_result = subprocess.run(
        [
            'nvidia-smi',
            '--query-compute-apps=gpu_uuid,pid,process_name,used_memory',
            '--format=csv,noheader,nounits'
        ],
        stdout=subprocess.PIPE,
        text=True,
        check=False,
    )

    if proc_result.stdout.strip():
        print("GPU Processes (with users):")
        for line in proc_result.stdout.strip().splitlines():
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 4:
                print(line)
                continue
            gpu_uuid, pid, proc_name, used_mem = parts[:4]
            user_result = subprocess.run(
                ['ps', '-o', 'user=', '-p', pid],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            user = user_result.stdout.strip() if user_result.returncode == 0 and user_result.stdout.strip() else 'unknown'
            print(f"{gpu_uuid}, PID={pid}, USER={user}, PROC={proc_name}, MEM={used_mem} MiB")
    else:
        print("GPU Processes (with users): none")
except Exception as e:
    print(f"Could not get GPU usage: {e}")
print(f"Best CUDA device ID: {best_dev_id}")