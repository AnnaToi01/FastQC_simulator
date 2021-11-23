from PIL import Image
# pillow required

image1 = Image.open(r'overrepresented.png')
image2 = Image.open(r'doublicates.png')

im1 = image1.convert('RGB')
im2 = image2.convert('RGB')

imagelist = [im2]

im1.save(r'result.pdf',save_all=True, append_images=imagelist)