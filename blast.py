import pandas as pd
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time
import glob

def run_blast_automation(input_csv, output_csv):
    # 1. CSV 파일 불러오기
    try:
        df = pd.read_csv(input_csv, header=None)
        df.columns = ['name', 'seq']
        print(f"총 {len(df)}개의 서열을 로드했습니다.")
    except Exception as e:
        print(f"파일을 읽는 중 오류 발생: {e}")
        return

    blast_results = []

    # 2. 각 행(row)에 대해 반복 수행
    for index, row in df.iterrows():
        seq_name = row['name']
        sequence = row['seq']
        
        print(f"[{index+1}/{len(df)}] '{seq_name}' BLAST 실행 중...", end="", flush=True)
        
        try:
            # 3. NCBI 서버에 BLAST 요청 (인터넷 연결 필수)
            # program="blastn" (DNA), database="nt" (nucleotide)
            # 단백질인 경우: program="blastp", database="nr" 사용
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
                    "e_value": hsp.expect,            # E-value (낮을수록 좋음)
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
            
            # 서버 과부하 방지를 위한 대기 (매너 모드)
            time.sleep(0.5) 

        except Exception as e:
            print(f" 에러 발생: {e}")
            blast_results.append({"query_name": seq_name, "error": str(e)})

    # 5. 결과를 DataFrame으로 변환 및 저장
    results_df = pd.DataFrame(blast_results)
    results_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    print(f"\n모든 작업 완료! 결과가 '{output_csv}'에 저장되었습니다.")

import argparse

# --- 실행 ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BLAST automation on a CSV file.")
    parser.add_argument("input_csv", help="Path to the input CSV file containing sequence data.")
    parser.add_argument("--output_csv", 
                        help="Path to the output CSV file for BLAST results. Defaults to input_csv_blast_results.csv",
                        default=None)

    args = parser.parse_args()

    input_file = args.input_csv
    output_file = args.output_csv

    if output_file is None:
        output_file = input_file.replace('.csv', '_blast_results.csv')

    print(f"Processing {input_file} -> {output_file}")
    run_blast_automation(input_file, output_file)