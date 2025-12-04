import pandas as pd
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO  # FASTA 파싱을 위해 추가
import time
import argparse
import os

def run_blast_fasta(input_fasta, output_csv):
    # 1. FASTA 파일 불러오기
    try:
        # SeqIO는 제너레이터이므로, 전체 개수 확인을 위해 리스트로 변환
        records = list(SeqIO.parse(input_fasta, "fasta"))
        print(f"총 {len(records)}개의 서열을 로드했습니다.")
    except Exception as e:
        print(f"파일을 읽는 중 오류 발생: {e}")
        return

    blast_results = []

    # 2. 각 서열(record)에 대해 반복 수행
    for index, record in enumerate(records):
        seq_name = record.id   # FASTA 헤더 (> 뒤의 이름)
        sequence = str(record.seq) # 서열 데이터 (문자열로 변환)
        
        print(f"[{index+1}/{len(records)}] '{seq_name}' BLAST 실행 중...", end="", flush=True)
        
        try:
            # 3. NCBI 서버에 BLAST 요청 (인터넷 연결 필수)
            # program="blastn" (DNA), database="nt" (nucleotide)
            result_handle = NCBIWWW.qblast("blastn", "nt", sequence, megablast=True)
            
            # 4. 결과 파싱 (XML 형식)
            blast_record = NCBIXML.read(result_handle)
            
            # 결과가 있는지 확인
            if blast_record.alignments:
                # 가장 첫 번째(가장 유사한) 결과만 가져오기
                alignment = blast_record.alignments[0]
                hsp = alignment.hsps[0] # High-scoring Segment Pair
                
                # 필요한 정보 추출
                result_data = {
                    "query_name": seq_name,
                    "subject_title": alignment.title, # 매칭된 서열의 이름
                    "length": alignment.length,       # 매칭된 서열 길이
                    "e_value": hsp.expect,            # E-value
                    "score": hsp.score,
                    "identity": f"{hsp.identities}/{hsp.align_length}", # 일치 개수
                    "identity_percentage": (hsp.identities / hsp.align_length) * 100
                }
                blast_results.append(result_data)
                print(" 완료!")
            else:
                print(" 결과 없음.")
                blast_results.append({"query_name": seq_name, "subject_title": "No hits found"})

            result_handle.close()
            
            # 서버 과부하 방지를 위한 대기
            time.sleep(1.0) # 안전하게 1초 대기

        except Exception as e:
            print(f" 에러 발생: {e}")
            blast_results.append({"query_name": seq_name, "error": str(e)})

    # 5. 결과를 DataFrame으로 변환 및 저장
    results_df = pd.DataFrame(blast_results)
    results_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    print(f"\n모든 작업 완료! 결과가 '{output_csv}'에 저장되었습니다.")

# --- 실행 ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BLAST automation on a FASTA file.")
    parser.add_argument("input_fasta", help="Path to the input FASTA file.")
    parser.add_argument("--output_csv", 
                        help="Path to the output CSV file for BLAST results.",
                        default=None)

    args = parser.parse_args()

    input_file = args.input_fasta
    output_file = args.output_csv

    # 출력 파일명이 지정되지 않았으면 입력 파일명 기반으로 생성
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_blast_results.csv"

    print(f"Processing {input_file} -> {output_file}")
    run_blast_fasta(input_file, output_file)