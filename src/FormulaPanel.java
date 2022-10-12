import java.awt.Font;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JPanel;

public class FormulaPanel extends JPanel implements MouseListener, KeyListener {

	private static final long serialVersionUID = 4631820944519889788L;
	static Images formulaBackground;

	// The label containing the String to be displayed
	private String label = "";
	private boolean singleElement = true, endBracket = false;

	// Calls the constructor, which adds keylistener and mouselistener
	FormulaPanel() throws IOException {
		// formulaBackground = new Images("FormulaPanel.png", 0, 0);
		this.addMouseListener(this);
		this.setVisible(true);
		this.addKeyListener(this);
		this.requestFocus();
	}

	// Displays the periodic table keyboard, the compound label, and its
	// subscripts
	@Override
	public void paintComponent(Graphics g) {

		g.drawImage(formulaBackground.img, formulaBackground.x, formulaBackground.y, null);
		Font n = new Font("SansSerif", Font.PLAIN, EquationPanel.fontFactor * 2);
		g.setFont(n);
		EquationPanel.compoundList.get(EquationPanel.compoundIndex).paintComponent(g, 510, 210);
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void mousePressed(MouseEvent e) {
		// TODO Auto-generated method stub

		int x = e.getX(), y = e.getY();
		int elementIndex = 0;

		// Checks to see which button the mouse clicked, and performs that
		// button's function
		if (x > 554 && x < 677 && y > 239 && y < 284 && singleElement) {
			Main.appState = 0;
			if (EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() > 0
					&& EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().get(
							EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() - 1)
							.getSubscript() == 0)
				EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList()
						.get(EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() - 1)
						.setSubscript(1);
			EquationPanel.compoundList.get(EquationPanel.compoundIndex).setMolarMass();
		} else if (EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() != 0 && x > 728
				&& x < 756 && y > 196 && y < 213) {
			singleElement = true;
			setLabel(
					getLabel()
							.substring(0,
									EquationPanel.compoundList.get(EquationPanel.compoundIndex)
											.getElementList().get(EquationPanel.compoundList
													.get(EquationPanel.compoundIndex).getElementList().size() - 1)
									.getIndex()));
			EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList()
					.remove(EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() - 1);
		} else if (x > 383 && x < 582 && y > 62 && y < 133) {
			if (EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() > 0
					&& EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().get(
							EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() - 1)
							.getSubscript() == 0)
				EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList()
						.get(EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() - 1)
						.setSubscript(1);
			EquationPanel.compoundList.get(EquationPanel.compoundIndex).setMolarMass();
			Main.appState = 2;
		} else if (y > 117 && y < 174) {
			if (x > 93 && x < 151) {
				elementIndex = 1;
			} else if (x > 1130 && x < 1188) {
				elementIndex = 2;
			}
		} else if (y > 178 && y < 235) {
			if (x > 93 && x < 151) {
				elementIndex = 3;
			} else if (x > 154 && x < 212) {
				elementIndex = 4;
			} else if (x > 826 && x < 884) {
				elementIndex = 5;
			} else if (x > 887 && x < 945) {
				elementIndex = 6;
			} else if (x > 948 && x < 1006) {
				elementIndex = 7;
			} else if (x > 1009 && x < 1067) {
				elementIndex = 8;
			} else if (x > 1070 && x < 1128) {
				elementIndex = 9;
			} else if (x > 1131 && x < 1189) {
				elementIndex = 10;
			}
		} else if (y > 238 && y < 296) {
			if (x > 93 && x < 151) {
				elementIndex = 11;
			} else if (x > 154 && x < 212) {
				elementIndex = 12;
			} else if (x > 826 && x < 884) {
				elementIndex = 13;
			} else if (x > 887 && x < 945) {
				elementIndex = 14;
			} else if (x > 948 && x < 1006) {
				elementIndex = 15;
			} else if (x > 1009 && x < 1067) {
				elementIndex = 16;
			} else if (x > 1070 && x < 1128) {
				elementIndex = 17;
			} else if (x > 1131 && x < 1189) {
				elementIndex = 18;
			}
		} else if (y > 299 && y < 540) {
			for (int i = 0; i < 18; i++) {
				if (x > i * 61 + 93 && x < i * 61 + 151) {
					elementIndex = i + 19 + 18 * ((y - 299) / 60);
					if (elementIndex > 56) {
						elementIndex += 14;
						if (elementIndex > 88) {
							elementIndex += 14;
						}
					}
					elementIndex = elementIndex == 71 || elementIndex == 103 || elementIndex > 109 ? 0 : elementIndex;
				}
			}
		} else if (y > 577 && y < 696) {
			for (int i = 0; i < 15; i++) {
				if (x > i * 60 + 215 && x < i * 60 + 273) {
					elementIndex = i + 57 + 32 * ((y - 577) / 60);
				}
			}
		}

		int sumElements = 0;

		for (int i = 0; i < EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size(); i++) {
			if (EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().get(i).isSingle()) {
				sumElements++;
			} else {
				sumElements += EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().get(i)
						.getElement().size();
			}
		}

		if (elementIndex != 0 && sumElements < 4) {
			Element element = null;
			try {
				element = new SingleElement(elementIndex, 0);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			ArrayList<Element> temp = EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList();
			element.setIndex(getLabel().length());
			if (temp.size() > 0 && temp.get(temp.size() - 1).getSubscript() < 1) {
				temp.get(temp.size() - 1).setSubscript(1);
			}

			if (singleElement) {
				element.setSingle(true);
				if (EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() < 4) {
					temp.add(element);
					setLabel(getLabel() + element.getAbbr() + " ");
				}
			} else if (temp.get(temp.size() - 1).getElement().size() < 4) {
				if (temp.get(temp.size() - 1).getElement().size() > 0 && temp.get(temp.size() - 1).getElement()
						.get(temp.get(temp.size() - 1).getElement().size() - 1).getSubscript() < 1) {
					temp.get(temp.size() - 1).getElement().get(temp.get(temp.size() - 1).getElement().size() - 1)
							.setSubscript(1);
				}
				temp.get(temp.size() - 1).getElement().add(element);
				setLabel(getLabel() + element.getAbbr() + " ");
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

	@Override
	public void keyTyped(KeyEvent e) {
		// TODO Auto-generated method stub

	}

	@Override
	public void keyPressed(KeyEvent e) {

	}

	// Sets a subscript if the inputted key is a number, creates a new Elemental
	// Compound if the inputed key is a '(', closes the most recent Elemental
	// Compound if the key is the ')'
	@Override
	public void keyReleased(KeyEvent e) {
		// TODO Auto-generated method stub

		int sumElements = 0;

		for (int i = 0; i < EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size(); i++) {
			if (EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().get(i).isSingle()) {
				sumElements++;
			} else {
				sumElements += EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().get(i)
						.getElement().size();
			}
		}

		char key = e.getKeyChar();
		if (key > 48 && key <= 57
				&& EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().size() > 0) {
			Compound temp = EquationPanel.compoundList.get(EquationPanel.compoundIndex);
			if (endBracket) {
				temp.getElementList().get(temp.getElementList().size() - 1).setSubscript(new Integer(key - 49));
				endBracket = false;
			} else if (singleElement) {
				temp.getElementList().get(temp.getElementList().size() - 1).setSubscript(key - 48);
			} else {
				temp.getElementList().get(temp.getElementList().size() - 1).getElement()
						.get(temp.getElementList().get(temp.getElementList().size() - 1).getElement().size() - 1)
						.setSubscript(key - 48);
			}
		} else if (key == '(' && singleElement && sumElements < 4) {
			Element temp = null;
			try {
				temp = new ElementalCompound();
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

			ArrayList<Element> elementList = EquationPanel.compoundList.get(EquationPanel.compoundIndex)
					.getElementList();
			temp.setIndex(getLabel().length());
			if (elementList.size() > 0 && elementList.get(elementList.size() - 1).getSubscript() < 1) {
				elementList.get(elementList.size() - 1).setSubscript(1);
			}
			EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList().add(temp);
			singleElement = false;
			setLabel(getLabel() + "(");
		} else if (key == ')' && !singleElement) {
			singleElement = true;
			endBracket = true;
			setLabel(getLabel() + ") ");
			ArrayList<Element> temp = EquationPanel.compoundList.get(EquationPanel.compoundIndex).getElementList();
			((ElementalCompound) (temp.get(temp.size() - 1))).setCompleted(true);
			temp.get(temp.size() - 1).setSubscript(0);
			if (temp.get(temp.size() - 1).getElement().size() > 0 && temp.get(temp.size() - 1).getElement()
					.get(temp.get(temp.size() - 1).getElement().size() - 1).getSubscript() < 1) {
				temp.get(temp.size() - 1).getElement().get(temp.get(temp.size() - 1).getElement().size() - 1)
						.setSubscript(1);
			}
			if (temp.get(temp.size() - 1).getElement().size() == 0) {
				temp.remove(temp.size() - 1);
				setLabel(getLabel().substring(0, getLabel().length() - "() ".length()));
			}

		}
	}

	// Accessors and mutators for this method 
	public String getLabel() {
		return label;
	}

	public void setLabel(String label) {
		this.label = label;
	}
}
