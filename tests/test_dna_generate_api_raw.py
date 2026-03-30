"""
DNA `/generate` API 纯净探针：脱离智能体，直接 POST 并打印原始 JSON。
运行: python tests/test_dna_generate_api_raw.py

可通过环境变量 DNA_GENERATE_API_URL 覆盖端点。
"""
import json
import os

import requests

API_URL = (os.getenv("DNA_GENERATE_API_URL") or "").strip() or "http://8.130.86.61:38008/generate"

# 与技能种子「官方示例」一致（服务文档示例 DNA）
payload = {
    "sequence": (
        "GAATAGGAACAGCTCCGGTCTACAGCTCCCAGCGTGAGCGACGCAGAAGACGGTGATTTCTGCATTTCCATCTGAGGTACCGGGTTCATCTCACTAGGGAGTGCCAGACAGTGGGCGCAGGCCAGTGTGTGTGCGCACCGTGCGCGAGCCGAAGCAGGGCGAGGCATTGCCTCACCTGGGAAGCGCAAGGGGTCAGGGAGTTCCCTTTCCGAGTCAAAGAAAGGGGTGATGGACGCACCTGGAAAATCGGGTCACTCCCACCCGAATATTGCGCTTTTCAGACCGGCTTAAGAAACGGCGCACCACGAGACTATATCCCACACCTGGCTCAGAGGGTCCTACGCCCACGGAATC"
    ),
    "num_tokens": 100,
    "temperature": 0.7,
    "top_k": 3,
    "top_p": 1.0,
    "enable_logits": False,
    "enable_sampled_probs": True,
}

if __name__ == "__main__":
    print("🚀 POST", API_URL)
    try:
        r = requests.post(API_URL, json=payload, headers={"Content-Type": "application/json"}, timeout=300)
        print("HTTP", r.status_code)
        try:
            print(json.dumps(r.json(), indent=2, ensure_ascii=False))
        except Exception:
            print(r.text[:4000])
    except Exception as e:
        print("❌", e)
