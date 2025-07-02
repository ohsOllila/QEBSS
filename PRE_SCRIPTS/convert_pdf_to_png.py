from pdf2image import convert_from_path

# Specify the path to your PDF file
pdf_path = 'ICL2_PRE_all.pdf'

# Convert PDF to a list of images (one for each page)
images = convert_from_path(pdf_path, dpi=300)

# Save each page as a PNG file
for i, image in enumerate(images):
    image.save('ICL2_PRE_all.png', 'PNG')

print("PDF has been converted to PNG images.")
