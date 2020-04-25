from PIL import Image, ImageDraw


im = Image.new('RGB', size=(300, 400))

x_origin = 20
y_origin = 20

cellSide = 50


y = y_origin

for j in range(5):
	x = x_origin
	for i in range(5):
		ImageDraw.Draw(im).polygon([(x,y),(x+cellSide,y),(x+cellSide,y+cellSide),(x,y+cellSide)], outline='white', fill='black')
		ImageDraw.Draw(im).line([(x+1,y+2),(x+cellSide-1,y+2)], fill='red', width=3)
		ImageDraw.Draw(im).line([(x+2,y+1),(x+2,y+cellSide-1)], fill='red', width=3)
		ImageDraw.Draw(im).ellipse([(x+15,y+15),(x+cellSide-15,y+cellSide-15)], fill='green')
		x = x + cellSide
	y = y + cellSide

x = x_origin
y = y_origin

ImageDraw.Draw(im).line([(x+1,y+2),(x+cellSide-1,y+2)], fill='red', width=3)

im.save('squares.png')