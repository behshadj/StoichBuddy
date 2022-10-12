import java.awt.Font;
import java.awt.Graphics;
import java.util.ArrayList;

public class Console {

	// No constructor since there is no instantiation 
	
	// An ArrayList of error messages
	static ArrayList<String> messageList = new ArrayList<String>();

	// Displays the error messages in the console, if there are more than 6
	// messages, it will delete the least recent one
	static void paintComponent(Graphics g) {
		Font n = new Font("SansSerif", Font.PLAIN, 13);
		g.setFont(n);

		if (messageList.size() > 6) {
			messageList.remove(0);
		}

		for (int i = 0; i < messageList.size(); i++) {
			g.drawString("> " + messageList.get(i), 913, 584 + i * 17);
		}
	}
}
