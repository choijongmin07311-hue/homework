import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

# Check if a filename is provided as a command-line argument
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = '20251127.csv' # Default filename

try:
    df = pd.read_csv(filename)
    original_data = df['close'].values
    original_data = original_data[:200]
    print(f"데이터 로드 성공 (데이터 길이: {len(original_data)}")
except FileNotFoundError:
    print("파일을 찾을 수 없습니다.")
    original_data = np.random.normal(0, 1, 200)

h1 = np.array([-1, -1, 0, 1, 1])
h2 = np.array([0.2, 0.2, 0.2, 0.2, 0.2])

y1 = np.convolve(original_data, h1, mode='same')
y2 = np.convolve(original_data, h2, mode='same')

plt.figure(figsize=(12, 10))

plt.subplot(3, 1, 1)
plt.title("stock price(close)")
plt.plot(original_data, color='black')
plt.grid(True)

plt.subplot(3, 1, 2)
plt.title("change detection (h1)")
plt.plot(y1, color='green')
plt.grid(True)

plt.subplot(3, 1, 3)
plt.title("moving average (h2)")
plt.plot(y2, color='blue')
plt.grid(True)

plt.tight_layout()
plt.savefig('result.jpg')