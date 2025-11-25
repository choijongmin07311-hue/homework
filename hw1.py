import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mping
import sys

filename = 'brain_mri.jpg'

try:
    original_img = mpimg.imread(filename)
    if len(original_img) == 3
        original_img = np.dot(original_img[...,:3], [0.2989, 0.5870, 0.1140])
    
except FileNotFoundError:
        print(f"파일 '{filename}'을 찾을 수 없습니다.")
        sys.exit()

mean_val = np.mean(original_img)
median_val= np.median(original_img)
print(f"이미지 평균값(Mean): {mean_val:.2f}")
print(f"이미지 중간값(Median): {median_val:.2f}")

img_mean_filtered[img_mean_filtered < mean_val] = -1
img_median_filtered[img_median_filtered < median] = -1

plt.figure(figsize=(15,5))

plt.subplot(1,3,1)
plt.title("Original Image")
plt.imshow(original_img, cmap='gray')
plt.axis('off')

