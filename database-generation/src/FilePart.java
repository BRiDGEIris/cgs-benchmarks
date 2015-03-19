import org.apache.commons.lang3.StringUtils;

public class FilePart implements Runnable {
	private int id = -1;
	private DatabaseExtension databaseExtension;
	private double ratio;

	private String[][] gt;
	private String[][] originalLines;

	private String[] sampleValues;
	private int currentAc[], currentAcHet[][], currentAcHetMerge[], currentAcHom[];	//The currentAcHetMerge is an intermediary value for easy computation, not present in the vcf
	private int currentAltQuantity;
	private String currentAlt[];

	private int combinatorialForAlt[];
	private long rand;

	public FilePart(int id, DatabaseExtension databaseExtension, String[][] originalLines) {
		this.id = id;
		this.databaseExtension = databaseExtension;
		this.originalLines = originalLines;
		this.sampleValues = new String[Main.SAMPLES];

		ratio = (Main.SAMPLES*1.0)/Main.TOTAL_SAMPLES;

		//We generate the available gt to reduce concatenation time below
		this.gt = new String[7][7];
		for(int i=0, j; i <= 6; i++)
			for(j=0; j <= 6; j++)
				gt[i][j] = "\t"+String.valueOf(i)+"/"+String.valueOf(j)+":";

		this.combinatorialForAlt = new int[7];
		this.combinatorialForAlt[0] = 0;
		this.combinatorialForAlt[1] = 1;
		this.combinatorialForAlt[2] = 3;
		this.combinatorialForAlt[3] = 6;
		this.combinatorialForAlt[4] = 10;
		this.combinatorialForAlt[5] = 15;
		this.combinatorialForAlt[6] = 21;
	}

	@Override
	public void run() {

		//We read each line, complete them and save them if we have a variant
		String currentLine;

		for(int lineId=0; lineId < originalLines.length && !originalLines[lineId][0].isEmpty(); lineId++) {	
			currentLine = this.completeLine(originalLines[lineId]);

			if(!currentLine.isEmpty())
				databaseExtension.setDestinationLine(id, lineId, currentLine.concat("\r\n"));			
		}
	}	

	/**
	 * Complete the current line
	 * @see http://demeranville.com/battle-of-the-tokenizers-delimited-text-parser-performance/
	 * @return
	 */
	private String completeLine(String[] elements){
		String tmp[], infos[],currentInfo[], currentSubInfo[], newLine, tmpText, dp;
		double oldProbAcHet[],oldProbAcHom[];
		int i, j, k, index, previousIndex, anAdj, altQuantity;
		//Split: ~10 times slower than indexOf... If the following lines would not be optimized the way
		//they are, you can wait 10 times longer if you want to have a few samples

		// We take some information from the INFO field. We need to count the quantity of ALT
		currentAlt = elements[4].split(",");
		altQuantity = currentAlt.length;

		//We need AN_Adj as the variants where not computed on all 60000 samples each time
		previousIndex = elements[7].indexOf("AN_Adj=") + 7;
		index = elements[7].indexOf(";", previousIndex);
		anAdj = Integer.valueOf(elements[7].substring(previousIndex, index));
		if(anAdj == 0)
			return "";
		ratio = (Main.SAMPLES*1.0)/(anAdj/2);		

		//We need DP
		previousIndex = elements[7].indexOf("DP=") + 3;
		index = elements[7].indexOf(";", previousIndex);
		dp = elements[7].substring(previousIndex, index);		

		//We also need the AC_Het and AC_Hom information
		previousIndex = elements[7].indexOf("AC_Het=") + 7;
		index = elements[7].indexOf(";", previousIndex);
		tmpText = elements[7].substring(previousIndex, index);
		tmp = tmpText.split(",");
		oldProbAcHet = new double[tmp.length];
		for(i=0; i < tmp.length; i++)
			oldProbAcHet[i] = Double.parseDouble(tmp[i]) / (anAdj/2);

		previousIndex = elements[7].indexOf("AC_Hom=") + 7;
		index = elements[7].indexOf(";", previousIndex);
		tmpText = elements[7].substring(previousIndex, index);
		tmp = tmpText.split(",");
		oldProbAcHom = new double[tmp.length];
		for(i=0; i < tmp.length; i++)
			oldProbAcHom[i] = Double.parseDouble(tmp[i]) / (anAdj/2);

		//We construct the samples information. currentAc & acHet & acHom are saved in attributes
		String sampleText = this.sampleValues(altQuantity, oldProbAcHet, oldProbAcHom, dp);
		if(sampleText.isEmpty())
			return "";

		//We create the alt column
		StringBuffer altTotal = new StringBuffer(currentAlt[0]);
		for(i=1; i < currentAlt.length && i < currentAltQuantity; i++)
			altTotal = altTotal.append(",").append(currentAlt[i]);

		//We construct the first part of the line
		newLine = elements[0]+"\t"+elements[1]+"\t"+elements[2]+"\t"+elements[3]+"\t"+altTotal+"\t"+elements[5]+"\t"+elements[6]+"\t";

		//We de-construct the INFO field, modify it, then reconstruct it
		infos = elements[7].split(";");
		int newAn = Main.SAMPLES * 2 + 2 * (int) Math.random()*Main.SAMPLES;
		for(i=0; i < infos.length; i++) {
			currentInfo = infos[i].split("=");
			currentSubInfo = currentInfo[0].split("_");

			if("AC".equals(currentInfo[0])) {

				//We simply give the AC values we have generated
				currentInfo[1] = "";
				for(j=0; j < this.currentAltQuantity; j++) {					
					if(!currentInfo[1].isEmpty())
						currentInfo[1] = currentInfo[1].concat(",");
					currentInfo[1] = currentInfo[1].concat(String.valueOf(this.currentAc[j]));			
				}

			} else if("AC".equals(currentSubInfo[0])) {

				if("Adj".equals(currentSubInfo[1])) {
					currentInfo[1] = "";
					for(j=0; j < this.currentAltQuantity; j++) {
						if(!currentInfo[1].isEmpty())
							currentInfo[1] = currentInfo[1].concat(",");
						currentInfo[1] = currentInfo[1].concat(String.valueOf(this.currentAcHetMerge[j] + this.currentAcHom[j]*2));
					}						
				} else if("Het".equals(currentSubInfo[1])) {
					currentInfo[1] = "";
					for(j=0; j < this.currentAltQuantity; j++) {
						for(k=j+1; k < this.currentAltQuantity+1; k++) {
							if(!currentInfo[1].isEmpty())
								currentInfo[1] = currentInfo[1].concat(",");
							currentInfo[1] = currentInfo[1].concat(String.valueOf(this.currentAcHet[j][k]));
						}
					}
				} else if("Hom".equals(currentSubInfo[1])) {
					currentInfo[1] = "";
					for(j=0; j < this.currentAltQuantity; j++) {
						if(!currentInfo[1].isEmpty())
							currentInfo[1] = currentInfo[1].concat(",");
						currentInfo[1] = currentInfo[1].concat(String.valueOf(this.currentAcHom[j]));
					}
				} else {
					// We couldn't take into account the AC_xx for the computation, so we simply estimate them
					tmp = currentInfo[1].split(",");
					currentInfo[1] = "";
					for(j=0; j < this.currentAltQuantity; j++) {
						if(!currentInfo[1].isEmpty())
							currentInfo[1] = currentInfo[1].concat(",");
						currentInfo[1] = currentInfo[1].concat(String.valueOf((int)(Integer.valueOf(tmp[j]) * ratio)));
					}
				}			
			} else if("AN".equals(currentSubInfo[0]) || "Het".equals(currentSubInfo[0]) || "Hom".equals(currentSubInfo[0])) {

				if(currentSubInfo.length > 1 && "Adj".equals(currentSubInfo[1])) {
					currentInfo[1] = String.valueOf(Main.SAMPLES * 2);					
				} else if("AN".equals(currentInfo[0])) {// We live in a perfect world with almost perfect machines
					currentInfo[1] = String.valueOf(newAn);
				} else {
					//We have no information about Het_xx, Hom_xx, AN_xx so we simply estimate them
					tmp = currentInfo[1].split(",");
					currentInfo[1] = "";
					for(j=0; j < tmp.length && j < this.combinatorialForAlt[this.currentAltQuantity]; j++) {
						if(!currentInfo[1].isEmpty())
							currentInfo[1] = currentInfo[1].concat(",");
						currentInfo[1] = currentInfo[1].concat(String.valueOf((int)(Integer.valueOf(tmp[j]) * ratio)));
					}
				}
			} else if("AF".equals(currentInfo[0])) {
				//We generate a new AF from the AC we have
				currentInfo[1] = "";
				for(j=0; j < this.currentAltQuantity; j++) {
					if(!currentInfo[1].isEmpty())
						currentInfo[1] = currentInfo[1].concat(",");
					currentInfo[1] = currentInfo[1].concat(String.valueOf(((double)currentAc[j])/newAn));
				}
			}

			//We merge the information
			if(i > 0)
				newLine = newLine.concat(";");

			if(currentInfo.length == 1)
				newLine = newLine.concat(currentInfo[0]);
			else
				newLine = newLine.concat(currentInfo[0]).concat("=").concat(currentInfo[1]);
		}		

		return newLine.concat(sampleText);
	}

	/**
	 * Construct the sample values. For now we only base our probabilities on the general AC_Het and AC_Hom (and 
	 * obviously AF), we do not take into account the Het_xx and Hom_xx. We will approximate them in the completeLine
	 * method. Note that we can also generate them here, if we got the information, without using them in the loop
	 * by simply using some random after we got the generated samples.
	 * TODO si AC_Het et AC_Hom = 0 je considère qu'on a aucune info, mais ce qui est bizarre
	 * c'est qu'on a AC > 0 avec AC_Het et AC_Hom = 0 dans la db initiale...
	 */
	private String sampleValues(int altQuantity, double[] oldProbAcHet, double[] oldProbAcHom, String dp) {
		double probHet = 0.0, probAdj = 0.0;
		double prob;
		String gt="";

		this.currentAc = new int[altQuantity];
		this.currentAcHet = new int[Math.min(6, altQuantity)][Math.min(6, altQuantity) + 1];
		this.currentAcHom =  new int[altQuantity];
		this.currentAcHetMerge = new int[altQuantity];
		int currentTotalAc = 0;
		double[][] pAcHet = new double[Math.min(6, altQuantity)][Math.min(6, altQuantity) + 1];

		double[] pAcHom = new double[altQuantity];
		pAcHom[0] = oldProbAcHom[0];
		for(int i=1; i < pAcHom.length; i++)
			pAcHom[i] = pAcHom[i-1] + oldProbAcHom[i];

		//We match x/y -> [x][y] for AC_Het values. Not sure about a clever algo so we do that stupidly
		//We directly compute cumulative probabilities, not the real acHet
		if(altQuantity == 1){
			pAcHet[0][1] = oldProbAcHet[0];
		} else if(altQuantity == 2) {
			pAcHet[0][1] = oldProbAcHet[0];
			pAcHet[0][2] = oldProbAcHet[1] + pAcHet[0][1];
			pAcHet[1][2] = oldProbAcHet[2] + pAcHet[0][2];
			//System.out.println("We give "+String.valueOf(oldProbAcHet[0]));
		} else if(altQuantity == 3) {
			pAcHet[0][1] = oldProbAcHet[0];
			pAcHet[0][2] = oldProbAcHet[1] + pAcHet[0][1];
			pAcHet[0][3] = oldProbAcHet[2] + pAcHet[0][2];
			pAcHet[1][2] = oldProbAcHet[3] + pAcHet[0][3];
			pAcHet[1][3] = oldProbAcHet[4] + pAcHet[1][2];
			pAcHet[2][3] = oldProbAcHet[5] + pAcHet[1][3];
		} else if(altQuantity == 4) {
			pAcHet[0][1] = oldProbAcHet[0];
			pAcHet[0][2] = oldProbAcHet[1] + pAcHet[0][1];
			pAcHet[0][3] = oldProbAcHet[2] + pAcHet[0][2];
			pAcHet[0][4] = oldProbAcHet[3] + pAcHet[0][3];
			pAcHet[1][2] = oldProbAcHet[4] + pAcHet[0][4];
			pAcHet[1][3] = oldProbAcHet[5] + pAcHet[1][2];
			pAcHet[1][4] = oldProbAcHet[6] + pAcHet[1][3];
			pAcHet[2][3] = oldProbAcHet[7] + pAcHet[1][4];
			pAcHet[2][4] = oldProbAcHet[8] + pAcHet[2][3];
			pAcHet[3][4] = oldProbAcHet[9] + pAcHet[2][4];
		} else if(altQuantity == 5) {
			pAcHet[0][1] = oldProbAcHet[0];
			pAcHet[0][2] = oldProbAcHet[1] + pAcHet[0][1];
			pAcHet[0][3] = oldProbAcHet[2] + pAcHet[0][2];
			pAcHet[0][4] = oldProbAcHet[3] + pAcHet[0][3];
			pAcHet[0][5] = oldProbAcHet[4] + pAcHet[0][4];
			pAcHet[1][2] = oldProbAcHet[5] + pAcHet[0][5];
			pAcHet[1][3] = oldProbAcHet[6] + pAcHet[1][2];
			pAcHet[1][4] = oldProbAcHet[7] + pAcHet[1][3];
			pAcHet[1][5] = oldProbAcHet[8] + pAcHet[1][4];
			pAcHet[2][3] = oldProbAcHet[9] + pAcHet[1][5];
			pAcHet[2][4] = oldProbAcHet[10] + pAcHet[2][3];
			pAcHet[2][5] = oldProbAcHet[11] + pAcHet[2][4];
			pAcHet[3][4] = oldProbAcHet[12] + pAcHet[2][5];
			pAcHet[3][5] = oldProbAcHet[13] + pAcHet[3][4];
			pAcHet[4][5] = oldProbAcHet[14] + pAcHet[3][5];
		} else if(altQuantity == 6) {
			pAcHet[0][1] = oldProbAcHet[0];
			pAcHet[0][2] = oldProbAcHet[1] + pAcHet[0][1];
			pAcHet[0][3] = oldProbAcHet[2] + pAcHet[0][2];
			pAcHet[0][4] = oldProbAcHet[3] + pAcHet[0][3];
			pAcHet[0][5] = oldProbAcHet[4] + pAcHet[0][4];
			pAcHet[0][6] = oldProbAcHet[5] + pAcHet[0][5];
			pAcHet[1][2] = oldProbAcHet[6] + pAcHet[0][6];
			pAcHet[1][3] = oldProbAcHet[7] + pAcHet[1][2];
			pAcHet[1][4] = oldProbAcHet[8] + pAcHet[1][3];
			pAcHet[1][5] = oldProbAcHet[9] + pAcHet[1][4];
			pAcHet[1][6] = oldProbAcHet[10] + pAcHet[1][5];
			pAcHet[2][3] = oldProbAcHet[11] + pAcHet[1][6];
			pAcHet[2][4] = oldProbAcHet[12] + pAcHet[2][3];
			pAcHet[2][5] = oldProbAcHet[13] + pAcHet[2][4];
			pAcHet[2][6] = oldProbAcHet[14] + pAcHet[2][5];
			pAcHet[3][4] = oldProbAcHet[15] + pAcHet[2][6];
			pAcHet[3][5] = oldProbAcHet[16] + pAcHet[3][4];
			pAcHet[3][6] = oldProbAcHet[17] + pAcHet[3][5];
			pAcHet[4][5] = oldProbAcHet[18] + pAcHet[3][6];
			pAcHet[4][6] = oldProbAcHet[19] + pAcHet[4][5];
			pAcHet[5][6] = oldProbAcHet[20] + pAcHet[4][6];
		} else {
			return "";
		}
		//Possible alternative to the ugly code above, but not sure if efficient:
		/*double lastValue = 0.0;
		for(int i=1, j=0, k=1; j < afValues.length; i++, k++) {
			pAcHet[j][i] = oldProbAcHet[k-1] + lastValue;
			lastValue = pAcHet[j][i];

			if(i == afValues.length) {
				j++;
				i=j;
			}				
		}*/



		//We compute the global probability to have a heterozygote and homozygote
		probHet = 0.0;
		for(int i=0; i < oldProbAcHet.length; i++) {
			if(i <= altQuantity)//We got 0/x or x/0, so 1 here = 1 sample.
				probHet += oldProbAcHet[i];
			else //We have x/y with x,y > 0. So 2 here = 1 sample.
				probHet += oldProbAcHet[i]/2.0;
		}
		probAdj = probHet;
		for(int i=0; i < oldProbAcHom.length; i++)
			probAdj += oldProbAcHom[i]; // We do not multiply oldAcHom by 2, as oldAcHom is related to the number of samples with no factor.

		/*
			We list each sample and generate gt & dp values. The instructions inside the loop
			will be executed 9 400 000 * 60 000 so be careful with what you put there. 
		 */

		int sampleLength = dp.length() + 1 + 3 + 2;
		StringBuffer text = new StringBuffer(sampleLength * Main.SAMPLES + 20);
		text = text.append("\tGT:DP");
		currentTotalAc = 0;
		this.XORShift64Init((long)Math.random()*Long.MAX_VALUE);
		for(int i=0,x=-1,y=-1; i < Main.SAMPLES; i++) {		
			prob = random01();//Math.random();

			if(prob <= probHet) {				
				//We have x/y or y/x with x > 0 or y > 0 or x,y >0. We have to find x & y.			
				outerloop:
				for(x=0; x < pAcHet.length; x++)
					for(y=x+1; y < pAcHet[x].length; y++)
						if(prob <= pAcHet[x][y])
							break outerloop;
				x = Math.min(x, 5);
				y = Math.min(y, altQuantity);

				//We decide the order x/y or y/x
				if(this.random01() < 0.5)
					gt = this.gt[x][y];
				else
					gt = this.gt[y][x];				
	
				this.currentAcHet[x][y] += 1;
	
				if(x > 0) {
					this.currentAc[x-1] += 1;
					this.currentAcHetMerge[x-1] += 1;
				}
	
				if(y > 0) {
					this.currentAc[y-1] += 1;//y is always > 0 but whatever
					this.currentAcHetMerge[y-1] += 1;
				}	
	
				currentTotalAc++;				
			} else if(prob <= probAdj) {
				//We have x/x with x > 0 and we have to find x.  
				prob -= probHet;				
				for(x=0; x < pAcHom.length; x++)
					if(prob <= pAcHom[x])
						break;
				x += 1;
				prob += probHet;


				gt = this.gt[x][x];
				this.currentAcHom[x-1] += 1;
				this.currentAc[x-1] += 1;

				currentTotalAc++;
			} else {
				gt = this.gt[0][0];
			}

			sampleValues[i] = new StringBuffer(gt).append(dp).toString();// -> 2 times faster than gt.concat(dp);	
		}	
		if(currentTotalAc == 0)
			return "";

		//We create the concatenated string		
		String str = StringUtils.join(sampleValues);

		//If we have more than one ALT possible, we need to check if each ALT is used, otherwise we need to change the informations...
		//Important information: we never have more than 6-7 possible ALT initially.
		if(altQuantity > 1) {
			int previousToReplace = -1;
			for(int altKey=1, j; altKey <= altQuantity; altKey++) {

				if(previousToReplace <= 0 && this.currentAcHom[altKey-1] == 0 && this.currentAcHetMerge[altKey-1] == 0) {
					previousToReplace = altKey;
				} else if(previousToReplace > 0 && !(this.currentAcHom[altKey-1] == 0 && this.currentAcHetMerge[altKey-1] == 0)) {
					//We only replace GT for the current ALT if the previous one was not used 
					str = StringUtils.replace(str, String.valueOf(altKey).concat("/"), String.valueOf(previousToReplace).concat("/"));
					str = StringUtils.replace(str, "/".concat(String.valueOf(altKey)), "/".concat(String.valueOf(previousToReplace)));

					this.currentAcHom[previousToReplace-1] = this.currentAcHom[altKey-1];
					this.currentAc[previousToReplace-1] = this.currentAc[altKey-1];
					this.currentAcHetMerge[previousToReplace-1] = this.currentAcHetMerge[altKey-1];
					this.currentAcHet[previousToReplace-1] = this.currentAcHet[altKey-1];
					for(j=0; j < altQuantity; j++)//A CHECKER avec tous les "-1" qui se rajoutent par-ci par-là je suis pas sûr
						this.currentAcHet[j][previousToReplace-1] = this.currentAcHet[j][altKey-1];		
					this.currentAlt[previousToReplace-1] = this.currentAlt[altKey-1];			

					this.currentAcHom[altKey-1] = 0;
					this.currentAcHetMerge[altKey-1] = 0;

					altKey = previousToReplace;//Don't worry.
					previousToReplace = -1;
				}
			}

			//We compute the quantity of alt used (way easier to do that here than above)
			this.currentAltQuantity = 0;
			for(int altKey=1; altKey <= altQuantity; altKey++)
				if(this.currentAcHom[altKey-1] > 0 || this.currentAcHetMerge[altKey-1] > 0)
					this.currentAltQuantity++;

		} else {
			this.currentAltQuantity = 1;
		}

		//We can finally return the complete line
		return text.append(str).toString();
	}


	/**
	 * Initialize a random generator better and faster (> 10x) than Math.random()
	 * @see http://demesos.blogspot.be/2011/09/replacing-java-random-generator.html
	 */
	public void XORShift64Init(long seed) {
		rand = seed==0 ? 0xdeadbeef : seed;
	}

	/**
	 * Return a random value better and faster (> 10x) than Math.random()
	 * @see http://demesos.blogspot.be/2011/09/replacing-java-random-generator.html
	 * @see http://en.wikipedia.org/wiki/Xorshift
	 */
	public double random01() {
		rand ^= (rand << 21); 
		rand ^= (rand >>> 35);
		rand ^= (rand << 4);
		
		return (rand*1.0 / Long.MAX_VALUE)/2.0 + 0.5;
	}
}
