# -*- coding: utf-8 -*-
"""2_Integrated Gradients.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1pTm-t8yRB5tnmLaxYvBEQcuFLyfWycT3

Follows https://www.tensorflow.org/tutorials/interpretability/integrated_gradients#calculate_integrated_gradients
"""

# %load_ext autoreload
# %load_ext autotime
# %autoreload 2

import os
import time
import math
import random
import numpy as np
import pandas as pd
from tqdm.notebook import tqdm
from pprint import pprint

import tensorflow as tf 
import PIL
from PIL import Image, ImageOps
import matplotlib.pyplot as plt
import pathlib
import shutil
import gc
from sys import getsizeof

# Mount my google drive
from google.colab import drive # will need to use verification code here
drive.mount('/content/gdrive', force_remount = True)

# Import all the images into the local environment i guess?
# May need to convert them to grayscale & size 512

def get_available_image_paths(data_root):
  #Get all image paths from my Drive
  all_image_paths = [str(path) for path in data_root.glob('**/*.tiff')]
  print('We have {} images'.format(len(all_image_paths)))
  print(all_image_paths[:2]) #get a sense of where the images are
  return all_image_paths

DEFAULT_IMG_SIZE = 512

def load_and_preprocess_image(path, img_size = DEFAULT_IMG_SIZE):
  #Can't use tensorflow preprocessing because tifs
  img = plt.imread(path)
  img = tf.image.resize(img, [img_size, img_size])
  img /= 255.0  # normalize pixels to 0,1
  return img

def convert_all_to_jpgs(source_folder, target_folder=None, delete_existing_jpg_folder=True, source_ext=".tiff"):

    target_folder = source_folder + "_jpg"
    if delete_existing_jpg_folder and os.path.exists(target_folder):
        shutil.rmtree(target_folder)
    os.mkdir(target_folder)

    for root, dirs, files in os.walk(source_folder, topdown=False):
        for dirname in dirs:
            os.mkdir(f'{target_folder}/{dirname}')
            
    for root, dirs, files in os.walk(source_folder, topdown=False):
        for name in files:
            infile = os.path.join(root, name)
            if '.tiff' in infile:
                outfile = infile.replace(source_folder, target_folder).replace('.tiff', ".jpg")
                with open(outfile, 'w+') as f:
                    im = Image.open(infile)
                    im.thumbnail(im.size)
                    im.save(f, "JPEG", quality=100)

root_path = 'gdrive/My Drive/Danino_Lab/Patterning_Scans/cheW_umoD_Expts/images_transformed_512'
data_root = pathlib.Path(root_path)
target_folder = 'images_transformed_512'
if not os.path.exists(target_folder):
  os.mkdir(target_folder)

# Get list of folders in the root folder & replicate in the local folder
for root, dirs, files in os.walk(data_root, topdown=False):
        for dirname in dirs:
          if not os.path.exists(f'{target_folder}/{dirname}'):
              os.mkdir(f'{target_folder}/{dirname}')
              print(dirname)

image_paths = get_available_image_paths(data_root)
for i, image_path in enumerate(image_paths):
  try:
    im = Image.open(image_path)
    # im.thumbnail((512, 512))
    img = im.resize((512, 512))
    gray_image = ImageOps.grayscale(img)
    # processed = load_and_preprocess_image(image_path, img_size = 512)
    newpath = image_path.replace('gdrive/My Drive/Danino_Lab/Patterning_Scans/cheW_umoD_Expts/images_transformed_512/', 'images_transformed_512/')
    gray_image.save(newpath)
  except Exception as e:
    print(e)
    print(image_path)
    break

# test what is going wrong with image size later??
image_path = image_paths[0]
im = Image.open(image_path)
print(im.size)
img = im.resize((512, 512))
print(img.size)
im.thumbnail((512, 512))
print(im.size)
gray_image = ImageOps.grayscale(im)
print(gray_image.size)
plt.imshow(gray_image)

# from image_utils import convert_all_to_jpgs
# if not os.path.exists('images_transformed_512_jpg'):
#     convert_all_to_jpgs('images_transformed_512')
# elif os.path.exists('images_transformed_512_jpg'):
#     data_dir = pathlib.Path('images_transformed_512_jpg')
#     if len(list((data_dir).glob('*/*.jpg'))) == 0:
#         convert_all_to_jpgs('images_transformed_512')
convert_all_to_jpgs('images_transformed_512')
data_dir = pathlib.Path('images_transformed_512_jpg')
print('We now have {} locally saved jpgs'.format(len(list(data_dir.glob('*/*.jpg')))))

# Check if all of the images are the wrong size?
for im_path in list(data_dir.glob('*/*.jpg')):
  im = Image.open(im_path)
  print(im.size)

# CHANGE THIS TO CHANGE IMAGE
# CHANGE THIS TO CHANGE IMAGE
# CHANGE THIS TO CHANGE IMAGE

chosen_fname = str(list(data_dir.glob('*/*.jpg'))[20])
chosen_image = PIL.Image.open(chosen_fname).convert('RGB')
chosen_image_tensor = tf.keras.utils.img_to_array(
    chosen_image, data_format=None, dtype=None
)
print(chosen_fname)
print(chosen_image_tensor.shape)
chosen_image

# Use a pure white image as a baseline, since our images are mostly white
baseline = tf.fill(chosen_image_tensor.shape, 255.0)

plt.imshow(baseline.numpy().astype('uint8'))
plt.title("Baseline")
plt.axis('off')
plt.show()

# Define image interpolation and visualize interpolation between our baseline and chosen image
def interpolate_images(baseline, image, alphas):
    alphas_x = alphas[:, tf.newaxis, tf.newaxis, tf.newaxis]
    baseline_x = tf.expand_dims(baseline, axis=0)
    input_x = tf.expand_dims(image, axis=0)
    delta = input_x - baseline_x
    images = baseline_x +  alphas_x * delta
    return images

m_steps=10 #50
alphas = tf.linspace(start=0.0, stop=1.0, num=m_steps+1) # Generate m_steps intervals for integral_approximation() below.

interpolated_images = interpolate_images(
    baseline=baseline,
    image=chosen_image_tensor,
    alphas=alphas)

fig = plt.figure(figsize=(20, 20))

i = 0
for alpha, image in zip(alphas[0::5], interpolated_images[0::5]):
    i += 1
    plt.subplot(1, len(alphas[0::5]), i)
    plt.title(f'alpha: {alpha:.1f}')
    plt.imshow(image.numpy().astype('uint8'))
    plt.axis('off')

plt.tight_layout();

print(interpolated_images.shape)

# Load our previously trained model (architecture + weights) and display its summary
# model = tf.keras.models.load_model('saved_inception_pretrained/model')
if not 'model' in locals():
  model = tf.keras.models.load_model('gdrive/My Drive/Danino_Lab/Patterning_Scans/cheW_umoD_Expts/saved_inception_pretrained/model')
model.summary()

getsizeof(model)

def compute_gradients(images, target_class_idx):
    with tf.GradientTape() as tape:
        tape.watch(images)
        probs = model(images)[:, target_class_idx]
    return tape.gradient(probs, images)

# only doing this to get the target index for our chosen image
train_dataset = tf.keras.utils.image_dataset_from_directory(
    data_dir,
    validation_split=0.5,
    subset="training",
    label_mode='categorical',
    seed=123,
    color_mode='rgb',
    image_size=(512, 512),
    batch_size=32)

class_names = train_dataset.class_names
target_class_idx = class_names.index(chosen_fname.split('/')[1])

path_gradients = compute_gradients(
    images=interpolated_images,
    target_class_idx=target_class_idx)
print(path_gradients.shape)

# Figure out this target_class_idx
print(target_class_idx)
print(class_names)
print(chosen_fname.split('/')[1])

def integral_approximation(gradients):
  # riemann_trapezoidal
    grads = (gradients[:-1] + gradients[1:]) / tf.constant(2.0)
    integrated_gradients = tf.math.reduce_mean(grads, axis=0)
    return integrated_gradients
ig = integral_approximation(gradients=path_gradients)
print(ig.shape)

@tf.function
def integrated_gradients(baseline,
                         image,
                         target_class_idx,
                         m_steps=50,
                         batch_size=32):
    # 1. Generate alphas.
    alphas = tf.linspace(start=0.0, stop=1.0, num=m_steps+1)

    # Initialize TensorArray outside loop to collect gradients.    
    gradient_batches = tf.TensorArray(tf.float32, size=m_steps+1)

    # Iterate alphas range and batch computation for speed, memory efficiency, and scaling to larger m_steps.
    for alpha in tf.range(0, len(alphas), batch_size):
        from_ = alpha
        to = tf.minimum(from_ + batch_size, len(alphas))
        alpha_batch = alphas[from_:to]

        # 2. Generate interpolated inputs between baseline and input.
        interpolated_path_input_batch = interpolate_images(baseline=baseline,
                                                           image=image,
                                                           alphas=alpha_batch)

        # 3. Compute gradients between model outputs and interpolated inputs.
        gradient_batch = compute_gradients(images=interpolated_path_input_batch,
                                           target_class_idx=target_class_idx)

        # Write batch indices and gradients to extend TensorArray.
        gradient_batches = gradient_batches.scatter(tf.range(from_, to), gradient_batch)    

    # Stack path gradients together row-wise into single tensor.
    total_gradients = gradient_batches.stack()

    # 4. Integral approximation through averaging gradients.
    avg_gradients = integral_approximation(gradients=total_gradients)

    # 5. Scale integrated gradients with respect to input.
    integrated_gradients = (image - baseline) * avg_gradients

    return integrated_gradients

ig_attributions = integrated_gradients(baseline=baseline,
                                       image=chosen_image_tensor,
                                       target_class_idx=target_class_idx,
                                       m_steps=10)
print(ig_attributions.shape)

# ORIGINAL CODE

def plot_img_attributions(baseline,
                          image,
                          target_class_idx,
                          attributions=None,
                          m_steps=50,
                          cmap=None,
                          overlay_alpha=1.0):
    print(f'Prediction class: {class_names[target_class_idx]}')
    if attributions is None:
        attributions = integrated_gradients(baseline=baseline, image=image, target_class_idx=target_class_idx, m_steps=m_steps)

    # Sum of the attributions across color channels for visualization.
    # The attribution mask shape is a grayscale image with height and width
    # equal to the original image.
    attribution_mask = tf.reduce_sum(tf.math.abs(attributions), axis=-1)

    fig, axs = plt.subplots(nrows=2, ncols=2, squeeze=False, figsize=(32, 32))

    axs[0, 0].set_title('Baseline image')
    axs[0, 0].imshow(baseline.numpy().astype('uint8'))
    axs[0, 0].axis('off')

    axs[0, 1].set_title('Original image')
    axs[0, 1].imshow(image.astype('uint8'))
    axs[0, 1].axis('off')

    axs[1, 0].set_title('Attribution mask')
    axs[1, 0].imshow(attribution_mask, cmap=cmap)
    axs[1, 0].axis('off')

    axs[1, 1].set_title('Overlay')
    axs[1, 1].imshow(attribution_mask, cmap=cmap)
    axs[1, 1].imshow(image.astype('uint8'), alpha=overlay_alpha)
    axs[1, 1].axis('off')

    plt.tight_layout()
    extent_baseline = axs[0, 0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    extent_original = axs[0, 1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    extent_attribution = axs[1, 0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    extent_overlay = axs[1, 1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())

    # Pad the saved area by 10% in the x-direction and 20% in the y-direction
    # fig.savefig('ax2_figure_expanded.png', bbox_inches=extent.expanded(1.1, 1.2))
    fig.savefig('baseline.svg', format='svg', bbox_inches=extent_baseline.expanded(2.0, 2.0))
    fig.savefig('original.svg', format='svg', bbox_inches=extent_original.expanded(2.0, 2.0))
    fig.savefig('attribution.svg', format='svg', bbox_inches=extent_attribution.expanded(2.0, 2.0))
    fig.savefig('overlay.svg', format='svg', bbox_inches=extent_overlay.expanded(2.0, 2.0))    
    return fig

plot_img_attributions(image=chosen_image_tensor,
                      baseline=baseline,
                      target_class_idx=target_class_idx,
                      attributions=ig_attributions,
                      m_steps=100,
                      cmap=plt.cm.inferno,
                      overlay_alpha=0.4)

# MODIFIED AD 091521
def plot_img_attributions(baseline,
                          image,
                          target_class_idx,
                          attributions=None,
                          m_steps=50,
                          cmap=None,
                          overlay_alpha=1.0):
    print(f'Prediction class: {class_names[target_class_idx]}')
    if attributions is None:
        attributions = integrated_gradients(baseline=baseline, image=image, target_class_idx=target_class_idx, m_steps=m_steps)

    # Sum of the attributions across color channels for visualization.
    # The attribution mask shape is a grayscale image with height and width
    # equal to the original image.
    attribution_mask = tf.reduce_sum(tf.math.abs(attributions), axis=-1)

    fig, axs = plt.subplots(nrows=2, ncols=2, squeeze=False, figsize=(32, 32))

    axs[0, 0].set_title('Baseline image')
    axs[0, 0].imshow(baseline.numpy().astype('uint8'))
    axs[0, 0].axis('off')

    axs[0, 1].set_title('Original image')
    axs[0, 1].imshow(image.astype('uint8'))
    axs[0, 1].axis('off')

    axs[1, 0].set_title('Attribution mask')
    axs[1, 0].imshow(attribution_mask, cmap=cmap)
    axs[1, 0].axis('off')

    axs[1, 1].set_title('Overlay')
    axs[1, 1].imshow(image.astype('uint8'))
    axs[1, 1].imshow(attribution_mask, cmap=cmap, alpha=overlay_alpha)

    # axs[1, 1].imshow(attribution_mask, cmap=cmap)
    # axs[1, 1].imshow(image.astype('uint8'), alpha=overlay_alpha)
    axs[1, 1].axis('off')

    plt.tight_layout()
    extent_baseline = axs[0, 0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    extent_original = axs[0, 1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    extent_attribution = axs[1, 0].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    extent_overlay = axs[1, 1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())

    # # Pad the saved area by 10% in the x-direction and 20% in the y-direction
    # # fig.savefig('ax2_figure_expanded.png', bbox_inches=extent.expanded(1.1, 1.2))
    # fig.savefig('baseline.svg', format='svg', bbox_inches=extent_baseline.expanded(2.0, 2.0))
    # fig.savefig('original.svg', format='svg', bbox_inches=extent_original.expanded(2.0, 2.0))
    # fig.savefig('attribution.svg', format='svg', bbox_inches=extent_attribution.expanded(2.0, 2.0))
    # fig.savefig('overlay.svg', format='svg', bbox_inches=extent_overlay.expanded(2.0, 2.0))    
    return fig

plot_img_attributions(image=chosen_image_tensor,
                      baseline=baseline,
                      target_class_idx=target_class_idx,
                      attributions=ig_attributions,
                      m_steps=100,
                      cmap=plt.cm.magma,
                      overlay_alpha=0.6)

# Sum of the attributions across color channels for visualization.
    # The attribution mask shape is a grayscale image with height and width
    # equal to the original image.
attribution_mask = tf.reduce_sum(tf.math.abs(ig_attributions), axis=-1)
print(tf.reduce_min(attribution_mask))
attribution_mask2 = attribution_mask
# attribution_mask2[attribution_mask==0] = 1.0
zerosmask = attribution_mask==0
print(zerosmask.dtype)
# attribution_mask2[zerosmask] = np.ones_like(zerosmask)
attribution_mask2 = np.where(attribution_mask==0, 0.02*np.ones(attribution_mask.shape), attribution_mask)
plt.imshow(attribution_mask2)
plt.colorbar()

alpha_mask = np.where(attribution_mask==0, 0.0, 0.4)
plt.imshow(chosen_image_tensor.astype('uint8'))
plt.imshow(attribution_mask, cmap=plt.cm.viridis, alpha=alpha_mask)
print(np.max(attribution_mask))
# plt.colorbar()

alpha_mask = np.where(attribution_mask==0, 0.0, 1.0)
plt.imshow(attribution_mask, cmap=plt.cm.plasma, alpha = alpha_mask)
plt.imshow(chosen_image_tensor.astype('uint8'), alpha = 0.4)

# let's try a nan approach instead, so that we can have a colorbar
attribution_mask2 = np.where(attribution_mask==0, np.nan, attribution_mask)
current_cmap = plt.cm.gist_heat_r
current_cmap.set_bad(color = 'white')
plt.imshow(attribution_mask2, cmap = current_cmap)
plt.imshow(chosen_image_tensor.astype('uint8'), alpha = 0.4)
plt.colorbar()

plt.figure()
plt.imshow(chosen_image_tensor.astype('uint8'))
plt.imshow(attribution_mask2, cmap = current_cmap, alpha = 0.6)
plt.colorbar()

# Show a zoomed in version
plt.figure()
plt.imshow(chosen_image_tensor[175:325, 175:325].astype('uint8'))
plt.imshow(attribution_mask2[175:325, 175:325], cmap = current_cmap, alpha = 0.3)
plt.colorbar()
plt.axis('off')
plt.title('overlay')

# Define the new plotting function

def plot_img_attributions(baseline,
                          image,
                          target_class_idx,
                          image_name,
                          attributions=None,
                          m_steps=50,
                          cmap=None,
                          overlay_alpha=0.4):
    # print(f'Prediction class: {class_names[target_class_idx]}')
    if attributions is None:
        attributions = integrated_gradients(baseline=baseline, image=image, target_class_idx=target_class_idx, m_steps=m_steps)

    # Sum of the attributions across color channels for visualization.
    # The attribution mask shape is a grayscale image with height and width
    # equal to the original image.
    attribution_mask = tf.reduce_sum(tf.math.abs(attributions), axis=-1)
    
    attribution_mask2 = np.where(attribution_mask==0, np.nan, attribution_mask)
    cmap.set_bad(color = 'white')
    fig = plt.figure()
    plt.imshow(attribution_mask2, cmap = cmap)
    plt.imshow(chosen_image_tensor.astype('uint8'), alpha = overlay_alpha)
    plt.colorbar()
    plt.title(image_name)
    plt.axis('off')
    plt.tight_layout()

    # fig.savefig(image_name+"overlay.svg", format='svg')    
    return fig

plot_img_attributions(image=chosen_image_tensor,
                      baseline=baseline,
                      target_class_idx=target_class_idx,
                      image_name = 'test',
                      attributions=ig_attributions,
                      m_steps=10,
                      cmap=plt.cm.gist_heat_r,
                      overlay_alpha=0.4)

# Let's see how many m_steps we can do before out of memory error, and if it looks diff w diff m_steps
plot_img_attributions(image=chosen_image_tensor,
                      baseline=baseline,
                      target_class_idx=target_class_idx,
                      image_name = 'test',
                      attributions=None,
                      m_steps=50,
                      cmap=plt.cm.gist_heat_r,
                      overlay_alpha=0.4)

impath = (list(data_dir.glob('*/*.jpg'))[20])
impathstr = str(impath)
print(impathstr)
print(impathstr.split('/')[2])
imname = impathstr.split('/')[2]
print(imname.replace('.jpg', '.svg'))
imname = imname.replace('Copy of Copy of ', '')
print(imname)

"""# Putting it all together: 
Let's use this method to make overlays for all the images, title each one as the image name, & save each
"""

# Make directories
source_folder = 'images_transformed_512_jpg'
target_folder = 'overlay_svgs'
if not os.path.exists(target_folder):
  os.mkdir(target_folder)
for root, dirs, files in os.walk(source_folder, topdown=False):
        for dirname in dirs:
            # print(f'{target_folder}/{dirname}')
            os.mkdir(f'{target_folder}/{dirname}')

num_ims = len(list(data_dir.glob('*/*.jpg')))
print(num_ims)

# Iterate over the image tensors
count = 0
for im_path in list(data_dir.glob('*/*.jpg')):
  # print the path?
  chosen_fname = str(im_path)
  # print(chosen_fname)

  # Get the image
  chosen_image = PIL.Image.open(chosen_fname).convert('RGB')
  
  chosen_image_tensor = tf.keras.utils.img_to_array(
      chosen_image, data_format=None, dtype=None)
  # Get the target class index
  target_class_idx = class_names.index(chosen_fname.split('/')[1])
  # Get the image name
  im_name = chosen_fname.split('/')[2]
  im_name = im_name.replace('Copy of Copy of ', str(count))

  # Make the attributions?
  img_attributions = integrated_gradients(baseline=baseline, image=chosen_image_tensor, target_class_idx=target_class_idx, m_steps=50)

  # Make the overlay figure
  fig = plot_img_attributions(image=chosen_image_tensor,
                      baseline=baseline,
                      target_class_idx=target_class_idx,
                      image_name = im_name,
                      attributions=img_attributions,
                      m_steps=50,
                      cmap=plt.cm.gist_heat_r,
                      overlay_alpha=0.4)
  # Title the overlay figure
  fig.title = im_name

  # Save the figure
  newpath = chosen_fname.replace('.jpg', 'overlay.svg')
  newpath = newpath.replace('images_transformed_512_jpg/', 'overlay_svgs/')
  newpath = newpath.replace('Copy of Copy of ', str(count))
  fig.savefig(newpath, format='svg')  

  # Clean up
  del chosen_image_tensor
  del chosen_image
  del img_attributions
  print(newpath)
  print(count)
  plt.close(fig)
  gc.collect()
  count+=1

# Zip the folder

# Zip the directory
shutil.make_archive('all_overlay_svgs', 'zip', 'overlay_svgs')

# this is a problem, everything has saved in the wrong place, now let's get a list of images, iterate over each, and move it
svgdir = pathlib.Path('/content')
print(len(list(svgdir.glob('*.svg'))))
print(len(list(svgdir.glob('*/*.svg'))))

for impath in list(svgdir.glob('*.svg')):
  impath = str(impath)
  # print(impath)
  newpath = impath.replace('/content', '/content/overlay_svgs')
  # print(newpath)
  shutil.move(impath, newpath)