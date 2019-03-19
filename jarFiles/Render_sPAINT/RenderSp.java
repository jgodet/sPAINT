import java.awt.Color;
import java.awt.Font;
import java.awt.image.IndexColorModel;
import java.util.ArrayList;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.TextRoi;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Duplicator;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


	
class RenderSp{
		private ImageProcessor ip;
		private ImageProcessor iprt;
		private ImageProcessor ipl;
		private ImageProcessor iplmin;
		private ColorProcessor cp;
		private IndexColorModel icm;
		private int height;
		private int width;
		private int[] frame;
		private boolean[] valid;
		private boolean[] validLambda;
		private double[] bkg ;
		private double[] signal;
		private double[] xcoord;
		private double[] ycoord;
		private double[] xsd;
		private double[] ysd;
		private boolean[] use;
		private double[] noise;
		private double[] lambda;
		private double[] precX;
		//private double[] precX2;
		private double[] precY;
		private double precMax;
		private int scaleFactor = 8;
		private double adu2phot = 1./55.; // to be modified
		private double pixSize = 106.67; //pixel size  
		private int boxSize =8;
		private int firstFrame;
		private int lastFrame;
		private String titlegd;
		private double waveMin;
		private double waveMax;



		public RenderSp(ImagePlus img, ArrayList<LocalisationRender> loc, int scaleFactor, double pixSize, double aduToPhoton, double precMax, int firstFrame, int lastFrame, String titlegd, double waveMin, double waveMax){
			
			/* Format ResultsTable op
			 * @ Frame(0)	origX(1)	origY(2)	origValue(3)	Error(4)	Noise(5)
			 * @ Background(6)	Signal(7)	Angle(8)	X(9)	Y(10)	X SD(11)	Y SD(12)	Precision(13)
			 * Z0Z1(14) lamba_nm(15)
			 * 
			 * à partir du 20150824
			 * @ X(0) Y(1) X.SD(2) Y.SD(3) Background(4) Noise(5) Signal(6) Frame(7) Lambda(8) Use(9)
			 * 
			 * 
			 * Noise - Background - Signal in ADU  
			 * X - Y - XSD - YSD - in pixel (uncalibrated)
			 * Precision in nm
			 * 
			 * 			
			 */
			
			float[] f = new float[loc.size()];
			float[] n = new float[loc.size()];
			float[] b = new float[loc.size()];
			float[] s = new float[loc.size()];
			float[] a = new float[loc.size()];
			float[] x = new float[loc.size()];
			float[] y = new float[loc.size()];
			float[] xs = new float[loc.size()];
			float[] ys = new float[loc.size()];
			float[] p = new float[loc.size()];
			float[] z = new float[loc.size()];
			float[] l = new float[loc.size()];
			boolean[] u = new boolean[loc.size()];  
			
			
			
			for (int i=0; i<loc.size(); i++){
				f[i]  = loc.get(i).getFrame();
				n[i]  = loc.get(i).getNoise();
				b[i]  = loc.get(i).getBkg();
				s[i]  = loc.get(i).getSignal();
				a[i]  = loc.get(i).getAngle();
				x[i]  = loc.get(i).getX();
				y[i]  = loc.get(i).getY();
				xs[i] = loc.get(i).getXsd();
				ys[i] = loc.get(i).getYsd();
				p[i]  = loc.get(i).getPrecision();
				z[i]  = loc.get(i).getZ0z1();
				l[i]  = loc.get(i).getLambda();
				u[i]  = loc.get(i).getUse();
			}

			this.height = (int) (img.getHeight() * scaleFactor);
			this.width = (int) (img.getWidth() * scaleFactor);
			
			//this.op = op;
			this.scaleFactor = scaleFactor;
			this.pixSize = pixSize;
			this.precMax = precMax;
			this.adu2phot = 1/aduToPhoton; 

			
			this.xcoord = scaleSize(convertFloatsToDoubles(x));
			this.ycoord = scaleSize(convertFloatsToDoubles(y));
			//this.xsd = scaleSD(convertFloatsToDoubles(xs));
			//this.ysd =  scaleSD(convertFloatsToDoubles(ys));
			this.xsd = convertFloatsToDoubles(xs);
			this.ysd = convertFloatsToDoubles(ys);
			this.bkg = convertToPhotons(convertFloatsToDoubles(b));
			this.noise = convertToPhotons(convertFloatsToDoubles(n));
			this.signal = convertToPhotons(convertFloatsToDoubles(s));
			this.frame = convertFloatsToInt(f);
			this.lambda =  convertFloatsToDoubles(l);
			this.use = u;
			
			//this.angle = convertFloatsToDoubles(op.getColumn(8));
			//this.prec = convertFloatsToDoubles(op.getColumn(13));
			//this.z0z1 =  convertFloatsToDoubles(op.getColumn(14));

			this.ip = new FloatProcessor(width, height);
			this.iprt = new FloatProcessor(width, height);
			this.ipl = new FloatProcessor(width, height);
			this.precX = new double[this.signal.length];
			this.precY = new double[this.signal.length];
			this.valid = new boolean[this.signal.length];
			this.validLambda = new boolean[this.signal.length];
			this.firstFrame = firstFrame;
			this.lastFrame = lastFrame;
			this.titlegd = titlegd;
			this.waveMin = waveMin;
			this.waveMax = waveMax;
			
			this.precX = sdPointing(this.xsd, this.signal, this.bkg, this.pixSize);
			//this.precX2 = sdPointing2(this.xsd, this.signal, this.bkg, this.pixSize);
			this.precY = sdPointing(this.ysd, this.signal, this.bkg, this.pixSize);

			validloc(this.frame, this.precX, this.precY);
			validlocLambda(this.frame, this.precX, this.precY, this.waveMin, this.waveMax);
			//TODO check diff entre les deux formules 
			/*for (int i=0; i<precX.length; i++){
				if(use[i]){
					System.out.println("Prec. Morten ="+precX[i] +";"+precX2[i] + " vs "+ p[i]);
				}
			}*/
		}
		
		public void render(){
			for (int i=0; i< xcoord.length; i++){
				if(valid[i]){
				intGauss(xcoord[i], ycoord[i], precX[i], precY[i]);
				}
			}
		}
		
		public void renderSpectral(){
			for (int i=0; i< xcoord.length; i++){
				if(validLambda[i] & use[i]){
				intGaussLambda(xcoord[i], ycoord[i], precX[i], precY[i], lambda[i]);
				}
			}
			lambdaCalc();
		}
		
		public void intGauss(double x, double y, double dx, double dy){
			/*
			 * calculate pixel values from Gaussian peak integrated over boxSize * boxSize pixels
			 * add the pixels values to an image processor ip (16 bits) used as an accumulator
			 * @x  : xloc (in pix)
			 * @y  : yloc (in pix)
			 * @dx : sd fitted gaussian x-axis
			 * @dy : sd fitted gaussian y-axis
			 * adapted from ....
			 */
			int xu = (int) x;
			int yv = (int) y;
			dx = dx * (double)scaleFactor / pixSize;
			dy = dy * (double)scaleFactor / pixSize;
			
			for (int u = xu - boxSize; u <= xu + boxSize; u++){
				for (int v = yv - boxSize; v <= yv + boxSize; v++){
					if( u > 0 & u < width){
						if (v > 0 & v < height){
							double xval =  (erf((u+1-x)/dx) - erf((u-x)/dx));
							double yval =  (erf((v+1-y)/dy) - erf((v-y)/dy));
							double val = xval*yval;
							iprt.setf(u,v,(float) (val+iprt.getf(u,v)));
						}
					}

		
				}
				
			}
		}
		
		public void intGaussLambda(double x, double y, double dx, double dy, double lambda){
			int xu = (int) x;
			int yv = (int) y;
			dx = dx * (double) scaleFactor / pixSize;
			dy = dy * (double) scaleFactor / pixSize;
			
			for (int u = xu - boxSize; u <= xu + boxSize; u++){
				for (int v = yv - boxSize; v <= yv + boxSize; v++){
					if( u > 0 & u < width){
						if (v > 0 & v < height){
							double xval = (erf((u+1-x) / dx) - erf((u-x) / dx));
							double yval = (erf((v+1-y) / dy) - erf((v-y) / dy));
							double  val = xval * yval;
							double vall = xval * yval * lambda;
							
							this.ip.setf(u, v, (float) (val+this.ip.getf(u,v)));	
							this.ipl.setf(u, v, (float) (vall+this.ipl.getf(u,v)));	
						}
					}

				}
				
			}
		}
		
		public void lambdaCalc(){
			for (int u = 0; u < width; u++ ){
				for (int v = 0; v < height; v++){
					if(this.ip.getf(u,v) > 0){
						double n = (double) this.ipl.getf(u,v);
						double d = (double) this.ip.getf(u,v);
						double l = n / d;
						this.ipl.setf(u,v,(float) (l));
					}
				}
			}		
		}
		
		public void applyLUT(double min, double max){
			this.ipl.setColorModel(icm);
			this.ipl.setMinAndMax(min, max);	
		}
		
		
		public ImageProcessor getIp(){
			return(ip);
		}

		public ImageProcessor getIpl(){
			return(ipl);
		}
		
		public ColorProcessor getCp(){
			return(cp);
		}
		
		public IndexColorModel getCm(){	
			return(icm);
		}
		
		public double erf(double z){
			/*
			 * from http://introcs.cs.princeton.edu/java/21function/ErrorFunction.java
			 * fractional error less than 1.2 * 10 ^ -7
			 * from Chebyshev fitting formula for erf(z) from Numerical Recipes, 6.2
			 * use Horner's method
			 */
			double t = 1.0 / (1.0 + .5 *Math.abs(z));
			
			double ans = 1 - t * Math.exp(-z * z - 1.265512223
					+ t * (1.00002368
					+ t * (0.37409196
					+ t * (0.09678418
				    + t * (-0.18628806
				    + t * (0.27886807
				    + t * (-1.13520398
				    + t * (1.48851587
				    + t * (-0.82215223
				    + t * (0.17087277)))))))))		
							);
			
					if( z >= 0 ){
						return ans;
					}else{
						return (-ans);
					}
		}
		
		public  double[] convertFloatsToDoubles(float[] input)
		{
		    if (input == null)
		    {
		        return null; // Or throw an exception 
		    }
		    double[] output = new double[input.length];
		    for (int i = 0; i < input.length; i++)
		    {
		        output[i] = input[i];
		    }
		    return output;
		}

		public  int[] convertFloatsToInt(float[] input)
		{
		    if (input == null)
		    {
		        return null; // Or throw an exception 
		    }
		    int[] output = new int[input.length];
		    for (int i = 0; i < input.length; i++)
		    {
		        output[i] = (int) input[i];
		    }
		    return output;
		}
		
		public double[] scaleSize(double[] input)
		{
		    if (input == null)
		    {
		        return null; // Or throw an exception 
		    }
		    double[] output = new double[input.length];
		    for (int i = 0; i < input.length; i++)
		    {
		        output[i] = (double)scaleFactor * input[i];
		    }
		    return output;
		}
		
		public double[] scaleSD(double[] input){
			if (input == null)
			 {
			        return null; // Or throw an exception 
			 }
			 double[] output = new double[input.length];
			 for (int i = 0; i < input.length; i++)
			 {
			   output[i] = (double)pixSize * input[i];
			 }
			    return output;
		}
		
		public  double[] convertToPhotons(double[] input)
		{
		    if (input == null)
		    {
		        return null; // Or throw an exception 
		    }
		    double[] output = new double[input.length];
		    for (int i = 0; i < input.length; i++)
		    {
		        output[i] =  adu2phot * input[i];
		    }
		    return output;
		}
		

		
		public void jet(int n){ // code adapted from ImageJ macro - J.Mutterer  
			/**
			 * Returns an IndexColorModel similar to MATLAB's jet color map.
			 *
			 * @param n the number of gray value (8-bit scale) to be used 
			 *            
			 * @return The "Jet" LUT with the specified number of color
			 * @see <a href="https://list.nih.gov/cgi-bin/wa.exe?A2=IMAGEJ;c8cb4d8d.1306">alternative
			 *      version</a> by Jérome Mutterer
			 */
			byte[] tabReds = new byte[n];
			byte[] tabGreens = new byte[n];
			byte[] tabBlues = new byte[n];
			tabReds[0] = tabGreens[0] = tabBlues[0] = 1; //fond noir
			for ( int i = 1; i < n; i++ ) {
				double i4= (double) (4.*i) / (double) (n);
				tabReds[i] = (byte) (255*Math.min(Math.max(Math.min(i4-1.5,-i4+4.5),0.),1.));
				tabGreens[i] = (byte) (255*Math.min(Math.max(Math.min(i4-0.5,-i4+3.5),0.),1.));
				tabBlues[i] = (byte) (255*Math.min(Math.max(Math.min(i4+0.5,-i4+2.5),0.),1.));
			}
			IndexColorModel cmjet = new IndexColorModel(8, n, tabReds, tabGreens, tabBlues);
			this.icm =  cmjet;
		}
		
		
		/*void applyJetLut(double min, double max) {
				
			// convert min and Max in 0-255 scale
			double maxipl = ipl.getMax();
			double minipl = ipl.getMin();
			double minLut = (min - minipl) *255. / (maxipl - minipl);
			if(minLut < 0){
				minLut = 0;
				}
			double maxLut = (max - minipl) * 255. / (maxipl - minipl);
			if(maxLut>255){
				maxLut = 255;
				}
					
			this.iplmin = ipl.convertToByte(true);
		
			// Display LUT	
			iplmin.setColorModel(icm);
			iplmin.setMinAndMax( minLut, maxLut);
		
			this.icm = (IndexColorModel) iplmin.getCurrentColorModel(); //update IndexColorModel aprÃ¨s setMinAndMax

		}  corrected June 20. 2016 JG*/
		
		
		void applyJetLut(double min, double max) {
			
			// convert min and Max in 0-255 scale		
			this.iplmin = ipl.convertToByte(true);
			for (int u = 0; u< ipl.getHeight(); u++){
				for (int v = 0; v < ipl.getWidth(); v++){
					double pix = (double)ipl.getf(u,v);
					if(pix >0){
						pix = (float) (255. * (float)(pix - min) / (max - min));
						if(pix<0){pix = 0;}
						if(pix>255){pix = 255;}
						iplmin.set(u,v,(int)pix);
					}
				}
			}
		
			// Display LUT	
			iplmin.setColorModel(icm);
			iplmin.setMinAndMax( 0, 255);
		
			this.icm = (IndexColorModel) iplmin.getCurrentColorModel(); //update IndexColorModel aprÃ¨s setMinAndMax

		}
		
		
		public void indexToRgbPlot(){ 
			
			// adapted from Burger & Burge p251 Index_To_Rgb.java
			int w = ipl.getWidth();
			int h = ipl.getHeight();
			
			int mapSize = this.icm.getMapSize();
			
			byte[] Rmap = new byte[mapSize]; this.icm.getReds(Rmap);
			byte[] Gmap = new byte[mapSize]; this.icm.getGreens(Gmap);
			byte[] Bmap = new byte[mapSize]; this.icm.getBlues(Bmap);
			
			double maxip = ip.getMax();
			double minip = ip.getMin();
			
			//create 24-bit RGB image
			ColorProcessor cp = new ColorProcessor(w,h);
			int[] RGB = new int[3];	
			float[] hsb = new float[3];
			
			for (int v = 0; v < h; v++){
				for (int u = 0; u < w; u++){
					int idx = iplmin.getPixel(u, v);
					double intens = (double) ip.getf(u, v); 
					float sat = (float) Math.max(0, Math.min((intens - minip) / (maxip - minip),1)); 
					
					int r = 0xff & Rmap[idx]; //unsigned bytes 
					int g = 0xff & Gmap[idx];
					int b = 0xff & Bmap[idx];
					RGB[0] = (byte) r; 
					RGB[1] = (byte) g; 
					RGB[2] = (byte) b; 
					
					hsb = Color.RGBtoHSB(r,g,b, hsb);
					int rgb = Color.HSBtoRGB(hsb[0], hsb[1], sat);
					cp.putPixel(u, v, rgb);
				
				}
			}
			
			this.cp = cp;
		}
		
		public void showIt(){
			ImagePlus cimg = new ImagePlus (""+titlegd+" RGB IMAGE", cp);
			cimg.show();
			//IJ.run(cimg, "Enhance Contrast", "saturated=0.35");
			Calibration cal =cimg.getCalibration();
			cal.setUnit("nm");
			double updatedPixSize = pixSize / scaleFactor;
			cal.pixelWidth = cal.pixelHeight = updatedPixSize;
			cal.xOrigin = cal.yOrigin = 0.;
			
			ImagePlus limg = new ImagePlus (""+titlegd+" Localisations", ip);
			limg.show();
			//IJ.run(limg, "Enhance Contrast", "saturated=0.35");
			Calibration calloc =limg.getCalibration();
			calloc.setUnit("nm");
			calloc.pixelWidth = calloc.pixelHeight = updatedPixSize;
			calloc.xOrigin = calloc.yOrigin = 0.;
			
			
		}

		public void showItTot(){
			
			ImagePlus cimg = new ImagePlus (""+titlegd+" RGB IMAGE", cp);
			cimg.show();
			//IJ.run(cimg, "Enhance Contrast", "saturated=0.35");
			Calibration cal =cimg.getCalibration();
			cal.setUnit("nm");
			double updatedPixSize = pixSize / scaleFactor;
			cal.pixelWidth = cal.pixelHeight = updatedPixSize;
			cal.xOrigin = cal.yOrigin = 0.;
			
			ImagePlus limg = new ImagePlus (""+titlegd+" Localisations", ip);
			limg.show();
			//IJ.run(limg, "Enhance Contrast", "saturated=0.35");
			Calibration calloc =limg.getCalibration();
			calloc.setUnit("nm");
			calloc.pixelWidth = calloc.pixelHeight = updatedPixSize;
			calloc.xOrigin = calloc.yOrigin = 0.;
			
			ImagePlus totimg = new ImagePlus (""+titlegd+" All Localisations", iprt);
			totimg.show();
			//IJ.run(totimg, "Enhance Contrast", "saturated=0.35");
			Calibration calloctot =totimg.getCalibration();
			calloctot.setUnit("nm");
			calloctot.pixelWidth = calloctot.pixelHeight = updatedPixSize;
			calloctot.xOrigin = calloctot.yOrigin = 0.;
			
			ImagePlus cimg2 =  new Duplicator().run(cimg);
			cimg2.show();
			IJ.run(cimg2, "Enhance Contrast", "saturated=0.35");
			cimg2.setTitle("a");
			
			
			ImagePlus imp = new ImagePlus (""+titlegd+" Overlay IMAGE", iprt);
			imp.show();
			IJ.run("Add Image...", "image=a x=0 y=0 opacity=100 zero");
			Calibration callocOv =imp.getCalibration();
			callocOv.setUnit("nm");
			callocOv.pixelWidth = callocOv.pixelHeight = updatedPixSize;
			callocOv.xOrigin = callocOv.yOrigin = 0.;
			
			cimg2.close();
		}

		
		public boolean[] validloc(int[] frame, double[] precX, double[] precY){    
			boolean a = true;
			for (int i =0; i < frame.length ; i++){
				a = precX[i] < precMax & precY[i] < precMax & frame[i]>=firstFrame & frame[i] <= lastFrame;
				valid[i] = a;	
			}
			return valid;
		}
		
		public boolean[] validlocLambda(int[] frame, double[] precX, double[] precY, double waveMin, double WaveMax){    
			boolean a = true;
		
			for (int i = 0; i < frame.length ; i++){
				if(use[i]){
					a = precX[i] < precMax & precY[i] < precMax & frame[i]>=firstFrame & frame[i] <= lastFrame & lambda[i] >= waveMin & lambda[i] <= waveMax;
				}else{
					a = false;
				}
					validLambda[i] = a;	
			}
			return validLambda;
		}

		/*public double[]  sdPointing2(double[] sdx, double[] signal, double[] bkg, double pixSize){
			// formula from J. Biomed. Opt. 17(4), 049801 (Apr 06, 2012). doi:10.1117/1.JBO.17.4.049801
			//TODO formula to be checked car donne pas memes résultats que celle Mortensen
			double[] precision = new double[this.signal.length];
			for (int i = 0; i < sdx.length; i++){
					precision[i] =  Math.sqrt(  ((2. * sdx[i] * sdx[i] + pixSize * pixSize/12.)/ (.95*signal[i])) + ((16. * Math.PI * sdx[i] * sdx[i] * sdx[i] * sdx[i] * .95 * bkg[i] ) / (signal[i] * signal[i] * .95 *.95 * pixSize * pixSize)));
			//QY detection of the camera set at 0.95 (almost constant over 550 to 650 nm)
			}
			return(precision);
		}*/
		
		public double[]  sdPointing(double[] sdx, double[] signal, double[] bkg, double pixSize){
			// formula from Mortensen 2010  Nat Methods. 2010 May ; 7(5): 377–381. doi:10.1038/nmeth.1447.
			
			double[] precision = new double[this.signal.length];
			for (int i = 0; i < sdx.length; i++){
				   double sigmaASq = sdx[i] * sdx[i] + pixSize * pixSize / 12. ; //excellent approximation for a \leq sdx
				   precision[i] = Math.sqrt( 2. * sigmaASq / signal[i] * (16. / 9. + (8 * Math.PI * sigmaASq * bkg[i] * bkg[i]) / (signal[i] * pixSize * pixSize)));
			
				   //QY detection of the camera set at 0.95 (almost constant over 550 to 650 nm)
			}
			return(precision);
		}
		
		public void drawLegend(double min, double max){
			FloatProcessor iplegend = new FloatProcessor(40 , 200 );
			for(int u = 0; u<40; u++){
				for (int v = 0; v < 200 ; v++){
					float val = (float)( max - (max-min) * (double) v / 200.) ;
					iplegend.setf(u,v,val);
				}
			}
			
			iplegend.setColorModel(icm);
			iplegend.setMinAndMax(min, max);
			
			ImageProcessor  iplegend8 = iplegend.convertToByte(true);
			ImagePlus plegend = new ImagePlus("Legend", iplegend8); 
			plegend.show();			
			IJ.run("Canvas Size...", "width="+ 80 +" height="+ 220 +" position=Center-Left");
			

			Font font = new Font("SansSerif", Font.PLAIN, 8);
			String stup = "" + max + " nm";
			String stmin = "" + min + " nm";
			String stmid = "" + (max + min)/2 + " nm";
			TextRoi roiup = new TextRoi(44, 6, stup, font); 
			TextRoi roimid = new TextRoi(44, 106, stmid, font);
			TextRoi roilow = new TextRoi(44, 200, stmin, font);
			//roi.drawPixels(fond);
			roiup.setStrokeColor(Color.white);
			roimid.setStrokeColor(Color.white);
			roilow.setStrokeColor(Color.white);
			Overlay ov = new Overlay(); 
			ov.add(roiup);
			ov.add(roimid);
			ov.add(roilow);
			plegend.setOverlay(ov);
		}		
}  // fin de la classe Rendering




	