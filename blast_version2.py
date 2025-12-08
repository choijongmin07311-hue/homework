import pandas as pd
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import time
import argparse
import os
import sys

# load_sequences 함수는 변경 사항 없음
def load_sequences(input_file):
    file_ext = os.path.splitext(input_file)[1].lower()
    seq_data = []

    try:
        if file_ext == '.csv':
            print("recognized as csv file")
            df = pd.read_csv(input_file, header=None)
            for _, row in df.iterrows():
                seq_data.append({'id': str(row[0]), 'seq': str(row[1])})
        
        elif file_ext in ['.fasta', '.fa', '.fna', 'fas']:
            print("recognized as fasta file")
            for record in SeqIO.parse(input_file, "fasta"):
                seq_data.append({'id': record.id, 'seq': str(record.seq)})
        
        else:
            print(f"cannot recognize file: {file_ext}")
            print("(only .csv, .fasta, .fa, .fna, fas are available")
            sys.exit(1)

        print(f"total {len(seq_data)} sequences are uploaded")
        return seq_data

    except Exception as e:
        print(f"fatal error is recognized: {e}")
        sys.exit(1)

def run_blast(input_file, output_csv):
    sequences = load_sequences(input_file)
    blast_results = []
    
    # 결과 개수 설정 (여기서 10개로 설정)
    NUM_HITS = 10 

    for index, item in enumerate(sequences):
        seq_name = item['id']
        sequence = item['seq']
        
        print(f"[{index+1}/{len(sequences)}] '{seq_name}' BLAST 실행 중...", end="", flush=True)
        
        try:
            # 1. hitlist_size 파라미터 추가: 상위 10개 결과를 요청
            result_handle = NCBIWWW.qblast("blastn", "nt", sequence, megablast=True, hitlist_size=NUM_HITS)
            
            blast_record = NCBIXML.read(result_handle)
            
            if blast_record.alignments:
                # 2. 결과가 있다면 반복문을 통해 여러 개의 Alignment를 처리
                # enumerate를 사용하여 순위(rank)도 같이 저장하면 좋습니다.
                for i, alignment in enumerate(blast_record.alignments):
                    
                    # 각 Alignment 안에도 여러 HSP(High-scoring Segment Pair)가 있을 수 있으나,
                    # 보통 가장 점수가 높은 첫 번째 HSP(hsps[0])를 대표로 사용합니다.
                    hsp = alignment.hsps[0]
                    
                    result_data = {
                        "query_name": seq_name,
                        "rank": i + 1,  # 1위부터 10위까지 순위 표시
                        "subject_title": alignment.title,
                        "length": alignment.length,
                        "e_value": hsp.expect,
                        "score": hsp.score,
                        "identity": f"{hsp.identities}/{hsp.align_length}",
                        "identity_percentage": (hsp.identities / hsp.align_length) * 100
                    }
                    blast_results.append(result_data)
                
                print(f" finish! ({len(blast_record.alignments)} hits found)")
            else:
                print(" no result")
                blast_results.append({"query_name": seq_name, "subject_title": "No hits found"})

            result_handle.close()
            # 서버 부하 방지를 위해 대기 시간을 조금 늘리는 것이 안전합니다 (10개를 긁어오므로)
            time.sleep(1.0) 

        except Exception as e:
            print(f" error: {e}")
            blast_results.append({"query_name": seq_name, "error": str(e)})

    results_df = pd.DataFrame(blast_results)
    results_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
    print(f"\n done: '{output_csv}'")

# --- 실행부 ---
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Auto-detect CSV/FASTA and run BLAST.")
    parser.add_argument("input_file", help="Input file path (.csv or .fasta)")
    parser.add_argument("--output_csv", help="Output CSV path", default=None)

    args = parser.parse_args()
    
    input_path = args.input_file
    output_path = args.output_csv

    if output_path is None:
        base_name = os.path.splitext(input_path)[0]
        output_path = f"{base_name}_blast_results.csv"

    print(f"Processing: {input_path} -> {output_path}")
    run_blast(input_path, output_path)
