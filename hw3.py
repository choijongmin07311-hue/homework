import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
import sys

filename = 'my_voice.wav'
try:
    fs, data = wavfile.read(filename)
    if len(data.shape) == 2:
        data = data.mean(axis=1)
    print(f"파일 '{filename}' 로드 성공! (주파수: {fs}Hz)")

except FileNotFoundError:
    print(f"오류: 파일 '{filename}'을 찾을 수 없습니다.")
    sys.exit()

data=data.astype(float)
energy = np.sum(data ** 2)
power = energy / len(data)
print(f"energy: {energy:.2f}")
print(f"power: {power:.2f}")

fft_result = np.fft.fft(data)
fft_magnitude =np.abs(fft_result)

freqs = np.fft.fftfreq(len(data), 1/fs)
half_idx = len(data) // 2
freqs_half = freqs[:half_idx]
mag_half = fft_magnitude[:half_idx]

plt.figure(figsize=(12, 8))

plt.subplot(2, 1, 1)
plt.title("Time Domain: Voice Waveform")
limit = int(fs * 0.05)
time_axis = np.arange(len(data))
plt.plot(time_axis[:limit], data[:limit])
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid(True)

plt.subplot(2, 1, 2)
plt.title("Frequency Domain: FFT Spectrum")
plt.plot(freqs_half, mag_half, color='red')
plt.xlim(0, 4000)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Magnitude")
plt.grid(True)

plt.tight_layout()
plt.savefig('output_plot.png')
plt.show()