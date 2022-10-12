import java.awt.image.BufferedImage;

import java.io.IOException;

import javax.imageio.ImageIO;

public class Images {

	// Declare the BufferedImage, x and y int values, a boolean determining if
	// the mouse is hovering over an image, and the name of the .png of the main
	// BufferedImage and the one assigned when the mouse is hovering over the
	// image
	BufferedImage img;
	int x, y;
	String name;

	// Declare the image's main image, x, and y values
	Images(String n, int imageX, int imageY) throws IOException {
		this.x = imageX;
		this.y = imageY;
		this.name = n;
		setImage(name);
	}
	
	// Assign the BufferedImage
	void setImage(String n) throws IOException {
		java.net.URL fileURL = getClass().getResource(n);
		img = ImageIO.read(fileURL);
	}

	// Checks 
	public boolean pointOnImage(int mouseX, int mouseY) {
		if (mouseX > x && mouseX < x + img.getWidth() 
		&& mouseY > y && mouseY < y + img.getHeight()) {
			return true;
		} else {
			return false;	
		}
		
	}

}
