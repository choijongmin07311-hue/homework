import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mping
import sys

if len(sys.argv) < 2:
    print("Usage: python hw1.py <filename>")
    sys.exit(1)
filename = sys.argv[1]

try:
    original_img = mping.imread(filename)
    if original_img.ndim == 3:
        original_img = np.dot(original_img[...,:3], [0.2989, 0.5870, 0.1140])
    original_img = original_img.astype(np.float64)
    
except FileNotFoundError:
        print(f"파일 '{filename}'을 찾을 수 없습니다.")
        sys.exit()

mean_val = np.mean(original_img)
median_val= np.median(original_img)
print(f"이미지 평균값(Mean): {mean_val:.2f}")
print(f"이미지 중간값(Median): {median_val:.2f}")

img_mean_filtered = original_img.copy()
img_median_filtered = original_img.copy()

img_mean_filtered[img_mean_filtered < mean_val] = -1
img_median_filtered[img_median_filtered < median_val] = -1

plt.figure(figsize=(15,5))

plt.subplot(1,3,1)
plt.title("Original Image")
plt.imshow(original_img, cmap='gray')
plt.axis('off')

plt.subplot(1,3,2)
plt.title(f"Filtered by Mean ({mean_val:.1f})")
plt.imshow(img_mean_filtered, cmap='gray')
plt.axis('off')

plt.subplot(1,3,3)
plt.title(f"Filtered by Median ({median_val:.1f})")
plt.imshow(img_median_filtered, cmap='gray')
plt.axis('off')
plt.tight_layout()
plt.savefig(f"'{filename}'_output.jpg")