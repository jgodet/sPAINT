import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.ArrayList;

import gdsc.core.utils.UnicodeReader;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
//import ij.measure.ResultsTable;
import ij.process.ShortProcessor;



public class Render_sPAINT implements PlugIn{
	/**
	 * This method gets called by ImageJ / Fiji.
	 *
	 * @param arg can be specified in plugins.config
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */

	  @Override
	  public void run(String arg) {

//public class Render_{
	/**
	 * This method gets called by ImageJ / Fiji.
	 *
	 * @param arg can be specified in plugins.config
	 * @see ij.plugin.PlugIn#run(java.lang.String)
	 */

	//public static void main(String[] args){
	//	String [] ij_args = {"-Dplugins.dir = /Applications/Fiji.app/plugins",
	//			            "-Dmacros.dir  = /Applications/Fiji.app/macros"};
	//	ij.ImageJ.main(ij_args);
		

		/*String path = IJ.getFilePath("Choose a Spectral Localisation File");
		ResultsTable op = new ResultsTable();	
		try{
			op = ResultsTable.open(path);
		} catch (IOException e){
			e.printStackTrace();
		}
		op.show("Spectral Localisations Table");
		*/
		//open the specral/localisation file
		ArrayList<LocalisationRender> localisations = new ArrayList<LocalisationRender>();

		
		 OpenDialog od = new OpenDialog("Open specral/localisation file...");
		    String fileName = od.getFileName();
		    String directory = od.getDirectory();
			if (fileName==null){
			  return;
			}
		
		BufferedReader input = null;
		 try{
					FileInputStream fis = new FileInputStream(directory + fileName);
					input = new BufferedReader(new UnicodeReader(fis, null));
					String line;
					while ((line = input.readLine()) != null){
						if (line.length() == 0)
							continue;
						if (line.charAt(0) == '#')
							continue;
	               //"frame,origX,origY,error,noise,bkg,signal,angle,x,y,xsd,ysd,precision,z0z1,lambda\n";
						String[] fields = line.split("\t");
						if (fields.length >= 6){
							int frame = Integer.parseInt(fields[0]);
							int origX = Integer.parseInt(fields[1]);
							int origY = Integer.parseInt(fields[2]);
							float error = Float.parseFloat(fields[3]);
							float noise = Float.parseFloat(fields[4]);
							float bkg = Float.parseFloat(fields[5]); //background
							float signal = Float.parseFloat(fields[6]);
							float angle = Float.parseFloat(fields[7]);
							float x = Float.parseFloat(fields[8]);
							float y = Float.parseFloat(fields[9]);
							float xsd = Float.parseFloat(fields[10]);
							float ysd = Float.parseFloat(fields[11]);
							float precision = Float.parseFloat(fields[12]);
							float z0z1 = Float.parseFloat(fields[13]);
							float lambda = Float.parseFloat(fields[14]);
							boolean use = Boolean.parseBoolean(fields[15]);
	
							localisations.add(new LocalisationRender(frame, origX,origY,error,noise,bkg,signal,angle,x,y,xsd,ysd,precision,z0z1, lambda, use));
						}
					}
				}
		    catch (IOException e){}
			finally{
				try{if (input != null) input.close();}
				catch (IOException e){}
			}
	 
		String titlegd = "";

			  GenericDialog gd = new GenericDialog("sPAINT Rendering");
			  gd.addStringField("Title: ", titlegd);
		      gd.addNumericField("Scale factor: ", 8, 0);
		      gd.addNumericField("Pixel Size (in nm): ", 107, 1);
		      gd.addNumericField("ADU to Photon CF", 62.1, 1);
		      gd.addNumericField("Precision Max (in nm)", 100, 1);
		      gd.addNumericField("First Frame", 1, 0);
		      gd.addNumericField("Last Frame", 1000, 0);
		      gd.addNumericField("Image Width", 512, 0);
		      gd.addNumericField("Image Heigth", 512, 0);
		      gd.addNumericField("Wavelength min", 450, 1);
		      gd.addNumericField("Wavelength max", 555, 1);
		      gd.showDialog();
		      if(gd.wasCanceled()) {
		    	  IJ.error("PlugIn canceled");
		    	  return;
		      }
		     
		      titlegd = gd.getNextString();
		      int scaleFact = (int) gd.getNextNumber();
		      double pixSize = gd.getNextNumber();
		      double aduToPhoton = gd.getNextNumber();
		      double precMax = gd.getNextNumber();
		      int firstFrame = (int) gd.getNextNumber();
		      int lastFrame = (int) gd.getNextNumber();
		      int width = (int)gd.getNextNumber();
		      int height = (int)gd.getNextNumber();
		      double waveMin = gd.getNextNumber();
		      double waveMax = gd.getNextNumber();
		      
		      if(firstFrame <1) firstFrame =1;

			  //int width = 512;
			  //int height = 512;
			  ShortProcessor ip = new ShortProcessor(width, height);
			  String title = "Reconstructed image";
    		  ImagePlus img = new ImagePlus(title, ip);
		      
		      RenderSp a = new RenderSp(img, localisations, scaleFact, pixSize, aduToPhoton, precMax, firstFrame, lastFrame, titlegd, waveMin, waveMax);
		      	
			  a.renderSpectral();
		      a.jet(256); //2 to 256 (8 bits)
			  a.applyJetLut(waveMin, waveMax);
			  a.indexToRgbPlot();
			  a.showIt();

	}
}

