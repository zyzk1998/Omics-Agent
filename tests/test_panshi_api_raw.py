"""
磐石 BepiPred3 API 纯净探针：脱离智能体，直接 POST 并打印原始 JSON。
运行: python tests/test_panshi_api_raw.py
"""
import json

import requests

API_URL = "http://120.220.102.26:38089/predict"

payload = {
    "url_or_content": ">7lj4_B\nRSTTLLALLALVLLYVSGALVFRALEQPHEQQAQRELGEVREKFLRAHPCVSDQELGLLIKEVADALGGGADPETQSTSHSAWDLGSAFFFSGTIITTIGYGNVALRTDAGRLFCIFYALVGIPLFGDILLAGVGDRLGSSLRHGIGHIEAIFLKWHVPPELVRVLSAEMLFLLIGCLLFVLTPTFVFCYMEDWSKLEAIYFVIVTLTTVGFGDYVAGADPRQDSPAYQPLVWFWILLGLPAYFASVLTTIGNWLRVVS",
    "top_epitope_percentage_cutoff": "top_20",
    "use_sequential_smoothing": False,
    "prediction_mode": "vt_pred",
}

if __name__ == "__main__":
    print("🚀 正在向磐石 API 发送请求，请耐心等待...")
    try:
        response = requests.post(API_URL, json=payload, timeout=300)
        print(f"✅ HTTP 状态码: {response.status_code}")

        raw_json = response.json()
        print("📦 原始 JSON 响应:")
        print(json.dumps(raw_json, indent=2, ensure_ascii=False))

    except Exception as e:
        print(f"❌ 请求发生异常: {e}")
