import java.util.ArrayList;
import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;


// usage:
// rungroovy cut_analysis.groovy <mc_data_file.hipo> <rga_data_file.hipo> <energy>
// for use on RGA data

double en = Double.parseDouble(args[2]);
double enmax = en+0.1; //GeV
double thetamax = 40;  //degrees
double phimax = 180;   //degrees
double vzmax = 50;
double wmin = 1.8;
double wmax = 0;

if(en > 7){wmax = 4.5;}
else if(en > 4){wmax = 4;}
else {wmax = 2.5;}

int bin_num = 50;

HipoDataSource reader = new HipoDataSource();

// The MC reconstructed 1D histos
H1F Q2_bin_1_mc = new H1F("Q2_bin_1_mc", "Q2_bin_1_mc", bin_num, 0,13);
H1F Q2_bin_2_mc = new H1F("Q2_bin_2_mc", "Q2_bin_2_mc", bin_num, 0,13);
H1F Q2_bin_3_mc = new H1F("Q2_bin_3_mc", "Q2_bin_3_mc", bin_num, 0,13);
H1F Q2_bin_4_mc = new H1F("Q2_bin_4_mc", "Q2_bin_4_mc", bin_num, 0,13);
H1F Q2_bin_5_mc = new H1F("Q2_bin_5_mc", "Q2_bin_5_mc", bin_num, 0,13);
H1F Q2_bin_6_mc = new H1F("Q2_bin_6_mc", "Q2_bin_6_mc", bin_num, 0,13);

H1F W_bin_1_mc = new H1F("W_bin_1_mc", "W_bin_1_mc", bin_num, wmin, wmax+0.5);
H1F W_bin_2_mc = new H1F("W_bin_2_mc", "W_bin_2_mc", bin_num, wmin, wmax+0.5);
H1F W_bin_3_mc = new H1F("W_bin_3_mc", "W_bin_3_mc", bin_num, wmin, wmax+0.5);
H1F W_bin_4_mc = new H1F("W_bin_4_mc", "W_bin_4_mc", bin_num, wmin, wmax+0.5);
H1F W_bin_5_mc = new H1F("W_bin_5_mc", "W_bin_5_mc", bin_num, wmin, wmax+0.5);
H1F W_bin_6_mc = new H1F("W_bin_6_mc", "W_bin_6_mc", bin_num, wmin, wmax+0.5);

H1F theta_bin_1_mc = new H1F("theta_bin_1_mc", "theta_bin_1_mc", bin_num, 0, thetamax+5);
H1F theta_bin_2_mc = new H1F("theta_bin_2_mc", "theta_bin_2_mc", bin_num, 0, thetamax+5);
H1F theta_bin_3_mc = new H1F("theta_bin_3_mc", "theta_bin_3_mc", bin_num, 0, thetamax+5);
H1F theta_bin_4_mc = new H1F("theta_bin_4_mc", "theta_bin_4_mc", bin_num, 0, thetamax+5);
H1F theta_bin_5_mc = new H1F("theta_bin_5_mc", "theta_bin_5_mc", bin_num, 0, thetamax+5);
H1F theta_bin_6_mc = new H1F("theta_bin_6_mc", "theta_bin_6_mc", bin_num, 0, thetamax+5);

// The reconstructed 1D histos 
H1F Q2_bin_1 = new H1F("Q2_bin_1", "Q2_bin_1", bin_num, 0, 13);
Q2_bin_1.setTitleX("Q2 [GeV^{2}]");

H1F Q2_bin_2 = new H1F("Q2_bin_2", "Q2_bin_2", bin_num, 0, 13);
Q2_bin_2.setTitleX("Q2 [GeV^{2}]");

H1F Q2_bin_3 = new H1F("Q2_bin_3", "Q2_bin_3", bin_num, 0, 13);
Q2_bin_3.setTitleX("Q2 [GeV^{2}]");

H1F Q2_bin_4 = new H1F("Q2_bin_4", "Q2_bin_4", bin_num, 0, 13);
Q2_bin_4.setTitleX("Q2 [GeV^{2}]");

H1F Q2_bin_5 = new H1F("Q2_bin_5", "Q2_bin_5", bin_num, 0, 13);
Q2_bin_5.setTitleX("Q2 [GeV^{2}]");

H1F Q2_bin_6 = new H1F("Q2_bin_6", "Q2_bin_6", bin_num, 0, 13);
Q2_bin_6.setTitleX("Q2 [GeV^{2}]");

H1F W_bin_1 = new H1F("W_bin_1", "W_bin_1", bin_num, wmin, wmax+0.5);
W_bin_1.setTitleX("W [GeV]");

H1F W_bin_2 = new H1F("W_bin_2", "W_bin_2", bin_num, wmin, wmax+0.5);
W_bin_2.setTitleX("W [GeV]");

H1F W_bin_3 = new H1F("W_bin_3", "W_bin_3", bin_num, wmin, wmax+0.5);
W_bin_3.setTitleX("W [GeV]");

H1F W_bin_4 = new H1F("W_bin_4", "W_bin_4", bin_num, wmin, wmax+0.5);
W_bin_4.setTitleX("W [GeV]");

H1F W_bin_5 = new H1F("W_bin_5", "W_bin_5", bin_num, wmin, wmax+0.5);
W_bin_5.setTitleX("W [GeV]");

H1F W_bin_6 = new H1F("W_bin_6", "W_bin_6", bin_num, wmin, wmax+0.5);
W_bin_6.setTitleX("W [GeV]");

H1F theta_bin_1 = new H1F("theta_bin_1", "theta_bin_1", bin_num, 0,thetamax+5);
theta_bin_1.setTitleX("theta [deg]");

H1F theta_bin_2 = new H1F("theta_bin_2", "theta_bin_2", bin_num, 0,thetamax+5);
theta_bin_2.setTitleX("theta [deg]");

H1F theta_bin_3 = new H1F("theta_bin_3", "theta_bin_3", bin_num, 0,thetamax+5);
theta_bin_3.setTitleX("theta [deg]");

H1F theta_bin_4 = new H1F("theta_bin_4", "theta_bin_4", bin_num, 0,thetamax+5);
theta_bin_4.setTitleX("theta [deg]");

H1F theta_bin_5 = new H1F("theta_bin_5", "theta_bin_5", bin_num, 0,thetamax+5);
theta_bin_5.setTitleX("theta [deg]");

H1F theta_bin_6 = new H1F("theta_bin_6", "theta_bin_6", bin_num, 0,thetamax+5);
theta_bin_6.setTitleX("theta [deg]");



double e_mass = 0.000511;
double p_mass = 0.93827203;
Vector3 zero = new Vector3(0.0, 0.0, 0.0);
LorentzVector p_vec = new LorentzVector();
p_vec.setVectM(zero, p_mass);
LorentzVector e_vec = new LorentzVector(0.0, 0.0, en, en);


// create for loop for all files
// open cooked.files
// read in line by line
// for each line, open and run analysis
// close file

// for MC data
new File('.', args[0]).eachLine { line ->
    reader.open(line);
    
    byte q_mc = 0;
    float weight = 0;
    int pid_mc = 0;
    float px_mc = 0;
    float py_mc = 0;
    float pz_mc = 0;
    float mom_mc = 0;
    double phi_mc = 0;
    double theta_mc =  0;
    double Q2_mc = 0; 
    double W_mc =0;
    double E_prime_mc = 0;
    double xB_mc = 0;
    
    byte sector = 0;
    int cal_row = 0;
    
    while (reader.hasEvent()) {
        DataEvent event = reader.getNextEvent();
        
         // get MC reconstructed data
       if ( event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter") ) {
            DataBank bank_mc = event.getBank("REC::Particle");
            DataBank bank_evn_mc = event.getBank("MC::Event");
            DataBank bank_cal_mc = event.getBank("REC::Calorimeter");
            
            
            for (int k = 0; k < bank_mc.rows(); k++) {
                // get values
                px_mc = bank_mc.getFloat("px", k);
                py_mc = bank_mc.getFloat("py", k);
                pz_mc = bank_mc.getFloat("pz", k);
                q_mc = bank_mc.getByte("charge", k);
                
                //weight = bank_evn_mc.getFloat("weight", 0);
                weight=1;
                // calculate values
                mom_mc = (float) Math.sqrt(px_mc * px_mc + py_mc * py_mc + pz_mc * pz_mc);
                phi_mc = Math.atan2((double) py_mc,(double) px_mc);
                theta_mc = Math.acos((double) pz_mc/(double) mom_mc);
                
                Vector3 e_vec_3_mc = new Vector3(px_mc, py_mc, pz_mc); //3 vector e'
                LorentzVector e_vec_prime_mc = new LorentzVector(); //4 vector e'
                e_vec_prime_mc.setVectM(e_vec_3_mc, e_mass);
                
                LorentzVector q_vec_mc = new LorentzVector(); //4 vector q
                q_vec_mc.copy(e_vec); //e - e'
                q_vec_mc.sub(e_vec_prime_mc);
                Q2_mc = -q_vec_mc.mass2(); //-q^2
                
                
                LorentzVector w_vec_mc = new LorentzVector(); //4 vector used to calculate W
                w_vec_mc.copy(p_vec); //p + q
                w_vec_mc.add(q_vec_mc);
                W_mc = w_vec_mc.mass();
                
                E_prime_mc = e_vec_prime_mc.e();
                xB_mc = Q2_mc/(2.0*p_mass*(en-E_prime_mc));
                
                //Q2_mc = 4*en*E_prime_mc*pow((sin(theta_mc/2.0)),2);
                
                theta_mc *= 180/Math.PI;
                phi_mc *= 180/Math.PI;
                
                if(q_mc != -1) continue;
                
                // cuts
                if (theta_mc < 5 || theta_mc > 40 || W_mc < 2 || Q2_mc < 1 || E_prime_mc < 0.1*en) { 
                   // System.out.println(xB_mc); 
                    continue; 
                }
                
                
                cal_row = cal_cut_row(event, k);
                //System.out.println(j + " " + bank_cal.rows());
                if(cal_row != -1){
                    sector = bank_cal_mc.getByte("sector",cal_row);
                    
                    float x_cal = bank_cal_mc.getFloat("x",cal_row);
                    float y_cal = bank_cal_mc.getFloat("y",cal_row);
                    
                    float lu = bank_cal_mc.getFloat("lu",cal_row);
                    float lv = bank_cal_mc.getFloat("lv",cal_row);
                    float lw = bank_cal_mc.getFloat("lw",cal_row);
                    
                    if(lu < 350 && lu > 60 && lv < 370 && lw < 390){
                        // binning
                        if(xB_mc < 0.15){ // && Q2_mc < 2){
                            Q2_bin_1_mc.fill(Q2_mc);
                            Q2_bin_1_mc.setLineColor(3); 
                            W_bin_1_mc.fill(W_mc);
                            W_bin_1_mc.setLineColor(3); 
                            theta_bin_1_mc.fill(theta_mc);
                            theta_bin_1_mc.setLineColor(3);
                        }
                        if(xB_mc > 0.15 && xB_mc < 0.2){ // && Q2_mc > 2 && Q2_mc < 2.5){
                            Q2_bin_2_mc.fill(Q2_mc);
                            Q2_bin_2_mc.setLineColor(3);
                            W_bin_2_mc.fill(W_mc);
                            W_bin_2_mc.setLineColor(3);
                            theta_bin_2_mc.fill(theta_mc);
                            theta_bin_2_mc.setLineColor(3);
                        }
                        if(xB_mc > 0.2 && xB_mc < 0.27){ // && Q2_mc > 2.5 && Q2_mc < 3){
                            Q2_bin_3_mc.fill(Q2_mc);
                            Q2_bin_3_mc.setLineColor(3);
                            W_bin_3_mc.fill(W_mc);
                            W_bin_3_mc.setLineColor(3);
                            theta_bin_3_mc.fill(theta_mc);
                            theta_bin_3_mc.setLineColor(3);
                        }
                        if(xB_mc > 0.27 && xB_mc < 0.4){ // && Q2_mc > 3 && Q2_mc < 3.75){
                            Q2_bin_4_mc.fill(Q2_mc);
                            Q2_bin_4_mc.setLineColor(3);
                            W_bin_4_mc.fill(W_mc);
                            W_bin_4_mc.setLineColor(3);
                            theta_bin_4_mc.fill(theta_mc);
                            theta_bin_4_mc.setLineColor(3);
                        }
                        if(xB_mc > 0.4 && xB_mc < 0.6){ // && Q2_mc > 3.75 && Q2_mc < 4.5){
                            Q2_bin_5_mc.fill(Q2_mc);
                            Q2_bin_5_mc.setLineColor(3);
                            W_bin_5_mc.fill(W_mc);
                            W_bin_5_mc.setLineColor(3);
                            theta_bin_5_mc.fill(theta_mc);
                            theta_bin_5_mc.setLineColor(3);
                        }
                        if(xB_mc > 0.6){ // && Q2_mc > 4.5){
                            Q2_bin_6_mc.fill(Q2_mc);
                            Q2_bin_6_mc.setLineColor(3);
                            W_bin_6_mc.fill(W_mc);
                            W_bin_6_mc.setLineColor(3);
                            theta_bin_6_mc.fill(theta_mc);
                            theta_bin_6_mc.setLineColor(3);
                        }
                    }
                }
                
            } // end for loop
        }// end if 
    } // end event
} // end open MC data file

// normalization of MC histos
Q2_bin_1_mc.normalize(Q2_bin_1_mc.integral());
Q2_bin_2_mc.normalize(Q2_bin_2_mc.integral());
Q2_bin_3_mc.normalize(Q2_bin_3_mc.integral());
Q2_bin_4_mc.normalize(Q2_bin_4_mc.integral());
Q2_bin_5_mc.normalize(Q2_bin_5_mc.integral());
Q2_bin_6_mc.normalize(Q2_bin_6_mc.integral());
W_bin_1_mc.normalize(W_bin_1_mc.integral());
W_bin_2_mc.normalize(W_bin_2_mc.integral());
W_bin_3_mc.normalize(W_bin_3_mc.integral());
W_bin_4_mc.normalize(W_bin_4_mc.integral());
W_bin_5_mc.normalize(W_bin_5_mc.integral());
W_bin_6_mc.normalize(W_bin_6_mc.integral());
theta_bin_1_mc.normalize(theta_bin_1_mc.integral());
theta_bin_2_mc.normalize(theta_bin_2_mc.integral());
theta_bin_3_mc.normalize(theta_bin_3_mc.integral());
theta_bin_4_mc.normalize(theta_bin_4_mc.integral());
theta_bin_5_mc.normalize(theta_bin_5_mc.integral());
theta_bin_6_mc.normalize(theta_bin_6_mc.integral());

// open the RGA data file
new File('.', args[1]).eachLine { line ->
    reader.open(line);
//    reader.open(args[0]);
    
    double emax = 0;
    phimax = 0;
    thetamax = 0;
    vzmax = 0;
    int counter = 0;
    byte sector = 0;
    int cal_row = 0;
    int dc_row = 0;
    float Q2_bin_min = 0;
    float Q2_bin_max = 0;
    
    // generated data variables
    float weight = 0;
    
    // reconstructed data variables
    int pid = 0;
    byte q = 0;
    float px = 0;
    float py = 0;
    float pz = 0;
    float beta =0;
    float mom = 0;
    double phi = 0;
    double theta =  0;
    float vz = 0;
    double Q2 = 0; 
    double W =0;
    double E_prime = 0;
    double xB = 0;
                
                
    while (reader.hasEvent()) {
        DataEvent event = reader.getNextEvent();
        
        double tot_xsect = 0;
        
        // get reconstructed data
        if (event.hasBank("REC::Particle") && event.hasBank("REC::Calorimeter") && event.hasBank("REC::Traj")) {
            DataBank bank_rec = event.getBank("REC::Particle");
            DataBank bank_cal = event.getBank("REC::Calorimeter");
            DataBank bank_traj = event.getBank("REC::Traj");
            
            //System.out.println("number of bank_rec rows: " + bank_rec.rows() + ", # of gen bank rows: " + bank_gen.rows() );
                
            for (int k = 0; k < bank_rec.rows(); k++) {
    
                pid = bank_rec.getInt("pid", k);
                q = bank_rec.getByte("charge", k);
                px = bank_rec.getFloat("px", k);
                py = bank_rec.getFloat("py", k);
                pz = bank_rec.getFloat("pz", k);
                beta = bank_rec.getFloat("beta", k);
    
                mom = (float) Math.sqrt(px * px + py * py + pz * pz);
                phi = Math.atan2((double) py,(double) px);
                theta = Math.acos((double) pz/(double) mom);
                
                vz = bank_rec.getFloat("vz", k);
    
                // pick electrons
                //if (pid != 11) continue;
                if(q != -1) continue;
                
                Vector3 e_vec_3 = new Vector3(px, py, pz); //3 vector e'
                LorentzVector e_vec_prime = new LorentzVector(); //4 vector e'
                e_vec_prime.setVectM(e_vec_3, e_mass);
                
                LorentzVector q_vec = new LorentzVector(); //4 vector q
                q_vec.copy(e_vec); //e - e'
                q_vec.sub(e_vec_prime);
                Q2 = -q_vec.mass2(); //-q^2
                
                LorentzVector w_vec = new LorentzVector(); //4 vector used to calculate W
                w_vec.copy(p_vec); //p + q
                w_vec.add(q_vec);
                W = w_vec.mass();
                
                E_prime = e_vec_prime.e();
                xB = Q2/(2.0*p_mass*(en-E_prime));
                
                //Q2 = 4*en*E_prime*pow((sin(theta/2.0)),2);
                
                theta *= 180/Math.PI;
                phi *= 180/Math.PI;
                
                // Calculate max values of each param
                if(e_vec_prime.e()>emax){emax = e_vec_prime.e();} 
                if(theta > thetamax){thetamax = theta;}
                if(phi > phimax){phimax = phi;}
                if(vz > vzmax){vzmax = vz;}
                
                // begin cuts
                if (theta < 5 || theta > 40) {continue;}  
                if (W < 2) {continue;}
                if (Q2 < 1) {continue;}
                //if (E_prime < 0.2*en) {continue;}
                
                cal_row = cal_cut_row(event, k);
                //System.out.println(j + " " + bank_cal.rows());
                if(cal_row != -1){
                    sector = bank_cal.getByte("sector",cal_row);
                    
                    float x_cal = bank_cal.getFloat("x",cal_row);
                    float y_cal = bank_cal.getFloat("y",cal_row);
                    
                    float lu = bank_cal.getFloat("lu",cal_row);
                    float lv = bank_cal.getFloat("lv",cal_row);
                    float lw = bank_cal.getFloat("lw",cal_row);
                    
                    
                    if(lu < 350 && lu > 60 && lv < 370 && lw < 390){
                        
                        // binning
                        if(xB < 0.15){ // && Q2 < 2){
                            Q2_bin_1.fill(Q2);
                            W_bin_1.fill(W);
                            theta_bin_1.fill(theta);
                        }
                        if(xB > 0.15 && xB < 0.2){ // && Q2 > 2 && Q2 < 2.5){
                            Q2_bin_2.fill(Q2);
                            W_bin_2.fill(W);
                            theta_bin_2.fill(theta);
                        }
                        if(xB > 0.2 && xB < 0.27){ // && Q2 > 2.5 && Q2 < 3){
                            Q2_bin_3.fill(Q2);
                            W_bin_3.fill(W);
                            theta_bin_3.fill(theta);
                        }
                        if(xB > 0.27 && xB < 0.4){ // && Q2 > 3 && Q2 < 3.75){
                            Q2_bin_4.fill(Q2);
                            W_bin_4.fill(W);
                            theta_bin_4.fill(theta);
                        }
                        if(xB > 0.4 && xB < 0.6){ // && Q2 > 3.75 && Q2 < 4.5){
                            Q2_bin_5.fill(Q2);
                            W_bin_5.fill(W);
                            theta_bin_5.fill(theta);
                        }
                        if(xB > 0.6){ // && Q2 > 4.5){
                            Q2_bin_6.fill(Q2);
                            W_bin_6.fill(W);
                            theta_bin_6.fill(theta);
                        }
                    }
                }
                
                dc_row = dc_cut_row(event, k);
                if(dc_row != -1){
                    float x_dc = bank_traj.getFloat("x",dc_row);
                    float y_dc = bank_traj.getFloat("y",dc_row);
                    float z_dc = bank_traj.getFloat("z",dc_row);
                    
                    double pos = Math.sqrt(x_dc*x_dc + y_dc*y_dc + z_dc*z_dc);
                    double theta_dc = Math.acos((double) z_dc/ pos);
                    double phi_dc = Math.atan2((double) y_dc,(double) x_dc);
                    
                    theta_dc *= 180/Math.PI;
                    phi_dc *= 180/Math.PI;
                    
                }
            } // end for
        } // end if    
    } // end while
    
} // end open file

// normalization of RGA histograms
Q2_bin_1.normalize(Q2_bin_1.integral());
Q2_bin_2.normalize(Q2_bin_2.integral());
Q2_bin_3.normalize(Q2_bin_3.integral());
Q2_bin_4.normalize(Q2_bin_4.integral());
Q2_bin_5.normalize(Q2_bin_5.integral());
Q2_bin_6.normalize(Q2_bin_6.integral());
W_bin_1.normalize(W_bin_1.integral());
W_bin_2.normalize(W_bin_2.integral());
W_bin_3.normalize(W_bin_3.integral());
W_bin_4.normalize(W_bin_4.integral());
W_bin_5.normalize(W_bin_5.integral());
W_bin_6.normalize(W_bin_6.integral());
theta_bin_1.normalize(theta_bin_1.integral());
theta_bin_2.normalize(theta_bin_2.integral());
theta_bin_3.normalize(theta_bin_3.integral());
theta_bin_4.normalize(theta_bin_4.integral());
theta_bin_5.normalize(theta_bin_5.integral());
theta_bin_6.normalize(theta_bin_6.integral());

boolean dc_cut(float X, float Y, int S){
    boolean result= false;
    if( (S==3 || S==4 || S==5 || (Y>X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y<X*Math.tan(Math.PI*((S-1)/3.0+1.0/9))))
    && (S==1 || S==2 || S==6 || (Y<X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y>X*Math.tan(Math.PI*((S-1)/3.0+1.0/9)))) ) result= true;
    
    return result;
}

int cal_cut_row(DataEvent event, int row){
    DataBank bank_cal = event.getBank("REC::Calorimeter");
    
    int row_index = 0;
    int cal_row_match = -1;
    
    for(int j = 0; j < bank_cal.rows(); j++){
        row_index = bank_cal.getInt("pindex",j);
        
        if(row_index == row){
            cal_row_match = j;
            break;
        }
    }
    return cal_row_match;
}

int dc_cut_row(DataEvent event, int row){
    DataBank bank_traj = event.getBank("REC::Traj");
    int row_index = 0;
    int det_id = 0;
    int cal_row_match = -1;
    for(int j = 0; j < bank_traj.rows(); j++){
        row_index = bank_traj.getInt("pindex",j);
        det_id = bank_traj.getInt("detId",j);
        
        if(row_index == row && det_id == 6){
            cal_row_match = j;
            break;
        }
    }
    return cal_row_match;
}

TCanvas can_1d_a = new TCanvas("can", 1100, 600);
can_1d_a.setTitle("W binned");
can_1d_a.divide(3,2);
can_1d_a.cd(0);
can_1d_a.draw(W_bin_1);
can_1d_a.draw(W_bin_1_mc,"same");
can_1d_a.getPad().setLegend(true);
can_1d_a.getPad().setLegendPosition(20, 20);
can_1d_a.cd(1);
can_1d_a.draw(W_bin_2);
can_1d_a.draw(W_bin_2_mc,"same");
can_1d_a.getPad().setLegend(true);
can_1d_a.getPad().setLegendPosition(20, 20);
can_1d_a.cd(2);
can_1d_a.draw(W_bin_3);
can_1d_a.draw(W_bin_3_mc,"same");
can_1d_a.getPad().setLegend(true);
can_1d_a.getPad().setLegendPosition(20, 20);
can_1d_a.cd(3);
can_1d_a.draw(W_bin_4);
can_1d_a.draw(W_bin_4_mc,"same");
can_1d_a.getPad().setLegend(true);
can_1d_a.getPad().setLegendPosition(20, 20);
can_1d_a.cd(4);
can_1d_a.draw(W_bin_5);
can_1d_a.draw(W_bin_5_mc,"same");
can_1d_a.getPad().setLegend(true);
can_1d_a.getPad().setLegendPosition(20, 20);
can_1d_a.cd(5);
can_1d_a.draw(W_bin_6);
can_1d_a.draw(W_bin_6_mc,"same");
can_1d_a.getPad().setLegend(true);
can_1d_a.getPad().setLegendPosition(20, 20);
can_1d_a.save("figs/bins/W_binned.png");

TCanvas can_1d_b = new TCanvas("can", 1100, 600);
can_1d_b.setTitle("theta binned");
can_1d_b.divide(3,2);
can_1d_b.cd(0);
can_1d_b.draw(theta_bin_1);
can_1d_b.draw(theta_bin_1_mc,"same");
can_1d_b.getPad().setLegend(true);
can_1d_b.getPad().setLegendPosition(20, 20);
can_1d_b.cd(1);
can_1d_b.draw(theta_bin_2);
can_1d_b.draw(theta_bin_2_mc,"same");
can_1d_b.getPad().setLegend(true);
can_1d_b.getPad().setLegendPosition(20, 20);
can_1d_b.cd(2);
can_1d_b.draw(theta_bin_3);
can_1d_b.draw(theta_bin_3_mc,"same");
can_1d_b.getPad().setLegend(true);
can_1d_b.getPad().setLegendPosition(20, 20);
can_1d_b.cd(3);
can_1d_b.draw(theta_bin_4);
can_1d_b.draw(theta_bin_4_mc,"same");
can_1d_b.getPad().setLegend(true);
can_1d_b.getPad().setLegendPosition(20, 20);
can_1d_b.cd(4);
can_1d_b.draw(theta_bin_5);
can_1d_b.draw(theta_bin_5_mc,"same");
can_1d_b.getPad().setLegend(true);
can_1d_b.getPad().setLegendPosition(20, 20);
can_1d_b.cd(5);
can_1d_b.draw(theta_bin_6);
can_1d_b.draw(theta_bin_6_mc,"same");
can_1d_b.getPad().setLegend(true);
can_1d_b.getPad().setLegendPosition(20, 20);
can_1d_b.save("figs/bins/theta_binned.png");

TCanvas can_1d_c = new TCanvas("can", 1100, 600);
can_1d_c.setTitle("Q2 binned");
can_1d_c.divide(3,2);
can_1d_c.cd(0);
can_1d_c.draw(Q2_bin_1);
can_1d_c.draw(Q2_bin_1_mc,"same");
can_1d_c.getPad().setLegend(true);
can_1d_c.getPad().setLegendPosition(20, 20);
can_1d_c.cd(1);
can_1d_c.draw(Q2_bin_2);
can_1d_c.draw(Q2_bin_2_mc,"same");
can_1d_c.getPad().setLegend(true);
can_1d_c.getPad().setLegendPosition(20, 20);
can_1d_c.cd(2);
can_1d_c.draw(Q2_bin_3);
can_1d_c.draw(Q2_bin_3_mc,"same");
can_1d_c.getPad().setLegend(true);
can_1d_c.getPad().setLegendPosition(20, 20);
can_1d_c.cd(3);
can_1d_c.draw(Q2_bin_4);
can_1d_c.draw(Q2_bin_4_mc,"same");
can_1d_c.getPad().setLegend(true);
can_1d_c.getPad().setLegendPosition(20, 20);
can_1d_c.cd(4);
can_1d_c.draw(Q2_bin_5);
can_1d_c.draw(Q2_bin_5_mc,"same");
can_1d_c.getPad().setLegend(true);
can_1d_c.getPad().setLegendPosition(20, 20);
can_1d_c.cd(5);
can_1d_c.draw(Q2_bin_6);
can_1d_c.draw(Q2_bin_6_mc,"same");
can_1d_c.getPad().setLegend(true);
can_1d_c.getPad().setLegendPosition(20, 20);
can_1d_c.save("figs/bins/Q2_binned.png");