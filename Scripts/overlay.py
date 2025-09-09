from PIL import Image

# Open the background image (circos.png)
background = Image.open('circos.png').convert('RGBA')

# Open the overlay image (legend_trial.png)
overlay = Image.open('legend.png').convert('RGBA')

# Get dimensions
bg_width, bg_height = background.size
ov_width, ov_height = overlay.size

# Calculate position to center the overlay
pos_x = (bg_width - ov_width) // 2
pos_y = (bg_height - ov_height) // 2

# Create a copy of the background to avoid altering the original
result = background.copy()

# Paste the overlay in the center with its alpha channel as mask
result.paste(overlay, (pos_x, pos_y), overlay)

# Save the result
result.save('circos_with_legend.png', format='PNG')
print("Overlay complete! Output saved as circos_with_legend.png.")
