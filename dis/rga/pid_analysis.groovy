import org.jlab.clas.physics.LorentzVector;
import org.jlab.clas.physics.Vector3;
import org.jlab.groot.data.H1F;
import org.jlab.groot.data.H2F;
import org.jlab.groot.ui.TCanvas;
import org.jlab.io.base.DataBank;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataSource;

// usage:
// rungroovy pid_analysis.groovy <mc_data.hipo> <filename.hipo> <energy>
// this is for use with RGA data

double en = Double.parseDouble(args[2]);
double enmax = en+0.1; //GeV
double thetamax = 40;  //degrees
double phimax = 180;   //degrees
double vzmax = 50;
double wmin = 0;
double wmax = 0;
if(en > 7){wmax = 4.5;}
else if(en > 4){wmax = 4;}
else {wmax = 2.5;}

HipoDataSource reader = new HipoDataSource();

// 1D MC spectra
H1F h_eprime_mc = new H1F("h_eprime_mc", "E'_mc", 20,0,9);
H1F h_epcal_mc = new H1F("h_epcal_mc", "E_{PCAL}_mc", 50,0,2);
H1F h_ecin_mc = new H1F("h_ecin_mc", "E_{ECin}_mc", 50,0,1);
H1F h_ecout_mc = new H1F("h_ecout_mc", "E_{ECout}_mc", 50,0,0.2);
H1F h_ectot_mc = new H1F("h_ectot_mc", "E_{ECtot}_mc", 50,0,1.5);
H1F h_nphe_mc = new H1F("h_nphe_mc", "nphe_mc", 50,0,50);

// 1D spectra
H1F h_eprime = new H1F("h_eprime", "E'", 20,0,9);
h_eprime.setTitleX("E' [GeV]");

H1F h_epcal = new H1F("h_epcal", "E_{PCAL}", 50,0,2);
h_epcal.setTitleX("E_{pcal} [GeV]");

H1F h_ecin = new H1F("h_ecin", "E_{ECin}", 50,0,1);
h_ecin.setTitleX("E_{ECin} [GeV]");

H1F h_ecout = new H1F("h_ecout", "E_{ECout}", 50,0,0.2);
h_ecout.setTitleX("E_{ECout} [GeV]");

H1F h_ectot = new H1F("h_ectot", "E_{ECtot}", 50,0,1.5);
h_ectot.setTitleX("E_{ECtot} [GeV]");

H1F h_nphe = new H1F("h_nphe", "nphe", 50,0,50);
h_nphe.setTitleX("nphe");

// 2D Histos
H2F EC_vs_PCAL = new H2F("EC_vs_PCAL", "EC_{tot} vs E_{pcal}", 100, 0, 1.2, 100, 0, 1);
EC_vs_PCAL.setTitleX("E_{PCAL} [GeV]");
EC_vs_PCAL.setTitleY("EC_{in} + EC_{out} [GeV]");

H2F Etot_vs_p = new H2F("Etot_vs_p", "E_{tot}/p vs p", 100, 0, 8, 100, 0.15, 1);
Etot_vs_p.setTitleX("p [GeV]");
Etot_vs_p.setTitleY("E_tot/p");


H2F Etot_vs_p_1 = new H2F("Etot_vs_p_1", "E_{tot}/p vs p cut", 50, 0, 8, 50, 0.15, 1);
Etot_vs_p_1.setTitleX("p [GeV]");
Etot_vs_p_1.setTitleY("E_tot/p");

double e_mass = 0.000511;
double p_mass = 0.93827203;
Vector3 zero = new Vector3(0.0, 0.0, 0.0);
LorentzVector p_vec = new LorentzVector();
p_vec.setVectM(zero, p_mass);
LorentzVector e_vec = new LorentzVector(0.0, 0.0, en, en);


// for MC data
new File('.', args[0]).eachLine { line ->
    reader.open(line);
    
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
    double E_prime_mc = 0;
    double xB = 0;
    float nphe = 0.0;
            
    while (reader.hasEvent()) {
        DataEvent event_mc = reader.getNextEvent();
        
        if (event_mc.getBank("RECHB::Cherenkov")){
            DataBank bank_cher_mc = event_mc.getBank("RECHB::Cherenkov");
        
            for(int k = 0; k < bank_cher_mc.rows(); k++){
                nphe = bank_cher_mc.getFloat("nphe", k);
                h_nphe_mc.fill(nphe);
                h_nphe_mc.setLineColor(3);
            } 
            
        }
        
        // get reconstructed data
        if (event_mc.hasBank("RECHB::Particle") && event_mc.hasBank("RECHB::Calorimeter") && event_mc.hasBank("REC::Traj")) {
            DataBank bank_rec = event_mc.getBank("RECHB::Particle");
            DataBank bank_cal = event_mc.getBank("RECHB::Calorimeter");
            DataBank bank_traj = event_mc.getBank("REC::Traj");
            DataBank ecal_hits = event_mc.getBank("ECAL::clusters");
       
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
                
                theta *= 180/Math.PI;
                phi *= 180/Math.PI;
                vz = bank_rec.getFloat("vz", k);
    
                // pick electrons
                //if (pid != 11) continue;
                if(q != -1) continue;
                
                Vector3 e_vec_3 = new Vector3(px, py, pz); //3 vector e'
                LorentzVector e_vec_prime = new LorentzVector(); //4 vector e'
                e_vec_prime.setVectM(e_vec_3, e_mass);
                cal_row = cal_cut_row(event_mc, k);
                
                LorentzVector q_vec = new LorentzVector(); //4 vector q
                q_vec.copy(e_vec); //e - e'
                q_vec.sub(e_vec_prime);
                Q2 = -q_vec.mass2(); //-q^2
                
                LorentzVector w_vec = new LorentzVector(); //4 vector used to calculate W
                w_vec.copy(p_vec); //p + q
                w_vec.add(q_vec);
                W = w_vec.mass();
                
                E_prime_mc = e_vec_prime.e();
                xB = Q2/(2.0*p_mass*(en-E_prime_mc));
                
                // begin cuts
                if (theta < 5 || theta > 40) {continue;}  
                if (W < 2) {continue;}
                if (E_prime_mc < 0.1*en) {continue;}
                if (Q2 < 1) {continue;}
               
                // Calorimeter cuts
                if(cal_row != -1){
                    float e_pcal_mc = 0;
                    float ec_in_mc = 0;
                    float ec_out_mc = 0;  
                    
                    sector = bank_cal.getByte("sector",cal_row);
                    
                    float x_cal = bank_cal.getFloat("x",cal_row);
                    float y_cal = bank_cal.getFloat("y",cal_row);
                    
                    float lu = bank_cal.getFloat("lu",cal_row);
                    float lv = bank_cal.getFloat("lv",cal_row);
                    float lw = bank_cal.getFloat("lw",cal_row);
                    
                    if(lu < 350 && lu > 60 && lv < 370 && lw < 390){
                    
                        if( bank_cal.rows() >= 3){
                            for (int n = 0; n < bank_cal.rows(); n++){
                                byte layer = bank_cal.getByte("layer",n);
                                float energy = bank_cal.getFloat("energy",n);
                                
                                if(layer == 1) e_pcal_mc += energy;
                                if(layer == 4) ec_in_mc += energy;
                                if(layer == 7) ec_out_mc += energy;
                                
                            }
                            float ec_tot_mc = ec_in_mc + ec_out_mc;
                            
                            h_eprime_mc.fill(E_prime_mc);
                            h_eprime_mc.setLineColor(3); 
                            h_epcal_mc.fill(e_pcal_mc);
                            h_epcal_mc.setLineColor(3);
                            h_ecin_mc.fill(ec_in_mc);
                            h_ecin_mc.setLineColor(3);
                            h_ecout_mc.fill(ec_out_mc);
                            h_ecout_mc.setLineColor(3);
                            h_ectot_mc.fill(ec_tot_mc);
                            h_ectot_mc.setLineColor(3);
                            
                        }
                    }
                }
            } // end for
        } // end if    
    } // end while
} // end open file

// create for loop for all files
// open cooked.files
// read in line by line
// for each line, open and run analysis
// close file
new File('.', args[1]).eachLine { line ->
    reader.open(line);
    
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
    float nphe = 0.0;
            
    while (reader.hasEvent()) {
        DataEvent event = reader.getNextEvent();
        
        if (event.getBank("RECHB::Cherenkov")){
            DataBank bank_cher = event.getBank("RECHB::Cherenkov");
        
            for(int k = 0; k < bank_cher.rows(); k++){
                nphe = bank_cher.getFloat("nphe", k);
                h_nphe.fill(nphe);
            } 
            
        }
        
        // get reconstructed data
        if (event.hasBank("RECHB::Particle") && event.hasBank("RECHB::Calorimeter") && event.hasBank("REC::Traj")) {
            DataBank bank_rec = event.getBank("RECHB::Particle");
            DataBank bank_cal = event.getBank("RECHB::Calorimeter");
            DataBank bank_traj = event.getBank("REC::Traj");
            DataBank ecal_hits = event.getBank("ECAL::clusters");
       
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
                
                theta *= 180/Math.PI;
                phi *= 180/Math.PI;
                vz = bank_rec.getFloat("vz", k);
    
                // pick electrons
                //if (pid != 11) continue;
                if(q != -1) continue;
                
                Vector3 e_vec_3 = new Vector3(px, py, pz); //3 vector e'
                LorentzVector e_vec_prime = new LorentzVector(); //4 vector e'
                e_vec_prime.setVectM(e_vec_3, e_mass);
                cal_row = cal_cut_row(event, k);
                
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
                
                // begin cuts
                if (theta < 5 || theta > 40) {continue;}  
                if (W < 2) {continue;}
                if (E_prime < 0.1*en) {continue;}
                if (Q2 < 1) {continue;}
               
                // Calorimeter cuts
                if(cal_row != -1){
                    float e_pcal = 0;
                    float ec_in = 0;
                    float ec_out = 0;  
                    
                    sector = bank_cal.getByte("sector",cal_row);
                    
                    float x_cal = bank_cal.getFloat("x",cal_row);
                    float y_cal = bank_cal.getFloat("y",cal_row);
                    
                    float lu = bank_cal.getFloat("lu",cal_row);
                    float lv = bank_cal.getFloat("lv",cal_row);
                    float lw = bank_cal.getFloat("lw",cal_row);
                    
                    if(lu < 350 && lu > 60 && lv < 370 && lw < 390){
                    
                        if( bank_cal.rows() >= 3){
                            for (int n = 0; n < bank_cal.rows(); n++){
                                byte layer = bank_cal.getByte("layer",n);
                                float energy = bank_cal.getFloat("energy",n);
                                
                                if(layer == 1) e_pcal += energy;
                                if(layer == 4) ec_in += energy;
                                if(layer == 7) ec_out += energy;
                                
                            }
                            float ec_tot = ec_in + ec_out;
                            
                            
                            h_eprime.fill(E_prime);
                            h_epcal.fill(e_pcal);
                            h_ecin.fill(ec_in);
                            h_ecout.fill(ec_out);
                            h_ectot.fill(ec_tot);
                            
                            EC_vs_PCAL.fill(e_pcal, ec_tot);
                            Etot_vs_p.fill(mom, (ec_tot+e_pcal)/mom);
                            
                            if(e_pcal > 0.1){
                             Etot_vs_p_1.fill(mom, (ec_tot+e_pcal)/mom);
                             }
                            
                            
                        }
                    }
                }
            } // end for
        } // end if    
    } // end while
} // end open file

boolean dc_cut(float X, float Y, int S){
    boolean result= false;
    if( (S==3 || S==4 || S==5 || (Y>X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y<X*Math.tan(Math.PI*((S-1)/3.0+1.0/9))))
    && (S==1 || S==2 || S==6 || (Y<X*Math.tan(Math.PI*((S-1)/3.0-1.0/9)) && Y>X*Math.tan(Math.PI*((S-1)/3.0+1.0/9)))) ) result= true;
    
    return result;
}

int cal_cut_row(DataEvent event, int row){
    DataBank bank_cal = event.getBank("RECHB::Calorimeter");
    
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

TCanvas can = new TCanvas("can", 1100, 800);
can.setTitle("Energy spectra");
can.divide(3,2);
can.cd(0);
can.draw(h_eprime);
can.draw(h_eprime_mc,"same");
can.getPad().setLegend(true);
can.getPad().setLegendPosition(20, 20);
can.cd(1);
can.draw(h_epcal);
can.draw(h_epcal_mc,"same");
can.getPad().setLegend(true);
can.getPad().setLegendPosition(20, 20);
can.cd(2);
can.draw(h_ecin);
can.draw(h_ecin_mc,"same");
can.getPad().setLegend(true);
can.getPad().setLegendPosition(20, 20);
can.cd(3);
can.draw(h_ecout);
can.draw(h_ecout_mc,"same");
can.getPad().setLegend(true);
can.getPad().setLegendPosition(20, 20);
can.cd(4);
can.draw(h_ectot);
can.draw(h_ectot_mc,"same");
can.getPad().setLegend(true);
can.getPad().setLegendPosition(20, 20);
can.cd(5);
can.draw(h_nphe);
can.draw(h_nphe_mc,"same");
can.getPad().setLegend(true);
can.getPad().setLegendPosition(20, 20);
can.save("figs/pid/energy_spectra.png");


TCanvas can1 = new TCanvas("can1", 800, 600);
can1.setTitle("E_{tot} vs E_{PCAL}");
can1.draw(EC_vs_PCAL);
can1.save("figs/pid/ECtot_vs_Epcal.png");

TCanvas can2 = new TCanvas("can2", 800, 600);
can2.setTitle("E_{tot}/p vs p");
can2.draw(Etot_vs_p);
can2.save("figs/pid/Etot_vs_p.png");

TCanvas can3 = new TCanvas("can3", 800, 600);
can3.setTitle("E_{tot}/p vs p");
can3.draw(Etot_vs_p_1);
can3.save("figs/pid/Etot_vs_p_1.png");