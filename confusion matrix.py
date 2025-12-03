import random

# 0 또는 1을 100번 뽑아서 리스트로 생성
y_true = [random.randint(0, 1) for _ in range(100)]
print("길이:", len(y_true))
print("내용:", y_true)

y_pred = [random.randint(0,1) for _ in range(100)]
print("길이:", len(y_pred))
print("내용:", y_pred)

TP = 0  # True Positive: 실제 1, 예측 1
TN = 0  # True Negative: 실제 0, 예측 0
FP = 0  # False Positive: 실제 0, 예측 1
FN = 0  # False Negative: 실제 1, 예측 0

for i in range(len(y_true)):
    actual = y_true[i]      # 실제값
    predict = y_pred[i]     # 예측값

    # 조건문을 통해 각 경우의 수 카운트
    if actual == 1 and predict == 1:
        TP += 1
    elif actual == 0 and predict == 0:
        TN += 1
    elif actual == 0 and predict == 1:
        FP += 1
    elif actual == 1 and predict == 0:
        FN += 1

# 4. 결과 출력
print(f"TP(진양성): {TP}")
print(f"TN(진음성): {TN}")
print(f"FP(위양성): {FP}")
print(f"FN(위음성): {FN}")

print("\n--- Confusion Matrix (행렬 형태) ---")
# sklearn과 같은 순서: [[TN, FP], [FN, TP]]
print(f"[{TN:2d}, {FP:2d}]")
print(f"[{FN:2d}, {TP:2d}]")