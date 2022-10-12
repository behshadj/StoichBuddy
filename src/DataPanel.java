import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

import javax.swing.JPanel;
import javax.swing.JTextField;

public class DataPanel extends JPanel implements MouseListener, KeyListener {

	static Images dataBackground;

	// Textboxes for the user to input data
	JTextField volume = new JTextField(20);
	JTextField mass = new JTextField(20);
	JTextField density = new JTextField(20);
	JTextField moles = new JTextField(20);

	// Calls the constructor, which puts the textboxes in their correct
	// positions
	public DataPanel() {
		this.addMouseListener(this);
		volume.setBounds(298, 290, 130, 26);
		this.add(volume);
		volume.addKeyListener(this);
		volume.requestFocus();
		mass.setBounds(298, 535, 130, 26);
		this.add(mass);
		density.setBounds(790, 290, 130, 26);
		this.add(density);
		moles.setBounds(790, 535, 130, 26);
		this.add(moles);
		this.setLayout(null);
	}

	// Displays the background
	@Override
	public void paintComponent(Graphics g) {
		g.drawImage(dataBackground.img, dataBackground.x, dataBackground.y, null);

	}

	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyPressed(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub
		System.out.println("X: " + e.getX());
		System.out.println("Y: " + e.getY() + "\n");
		int x = e.getX(), y = e.getY();

		// return to FormulaPanel if the user pressed the "Chemical Formula"
		// button
		if (x > 651 && x < 867 && y > 66 && y < 140) {
			if (setData()) {
				Main.appState = 1;
			}
		}
		// return to EquationPanel if the user pressed the "Enter" putton
		else if (x > 555 && x < 676 && y > 374 && y < 417) {
			if (setData()) {
				Main.appState = 0;
			}
		}
	}

	@Override
	public void mouseReleased(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseEntered(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mouseExited(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	// Checks each data textbox (volume, mass, density, moles) to make sure
	// input is proper (i.e. numbers) and then passes their values to the
	// corresponding ones in the compound class. Returns false if an error
	// occurs, and true if all the data is successfully passed to the compound
	public boolean setData() {
		Compound compound = EquationPanel.compoundList.get(EquationPanel.compoundIndex);

		// Trys and catches are used because if an error occurs, the data in the
		// boxes simply resets

		try {
			if (volume.getText().indexOf('.') == -1 && volume.getText().length() > 0) {
				compound.setVolume((double) Integer.parseInt(volume.getText()));
			} else if (volume.getText().length() > 0) {
				compound.setVolume(Double.parseDouble(volume.getText()));
			}
		} catch (NumberFormatException n) {
			volume.setText("");
			return false;
		}

		try {
			if (mass.getText().indexOf('.') == -1 && mass.getText().length() > 0) {
				compound.setMass((double) Integer.parseInt(mass.getText()));
			} else if (mass.getText().length() > 0) {
				compound.setMass(Double.parseDouble(mass.getText()));
			}
		} catch (NumberFormatException n) {
			mass.setText("");
			return false;
		}

		try {
			if (density.getText().indexOf('.') == -1 && density.getText().length() > 0) {
				compound.setDensity((double) Integer.parseInt(density.getText()));
			} else if (density.getText().length() > 0) {
				compound.setDensity(Double.parseDouble(density.getText()));
			}
		} catch (NumberFormatException n) {
			density.setText("");
			return false;
		}

		try {
			if (moles.getText().indexOf('.') == -1 && moles.getText().length() > 0) {
				compound.setMoles((double) Integer.parseInt(moles.getText()));
			} else if (moles.getText().length() > 0) {
				compound.setMoles(Double.parseDouble(moles.getText()));
			}
		} catch (NumberFormatException n) {
			moles.setText("");
			return false;
		}
		return true;
	}
}
