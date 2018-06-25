#include "util.h"
#include "constants.h"
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cstdio>

#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
using namespace std;

double m_wd_array2[] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2};
static vector<double> m_wd_array(m_wd_array2, m_wd_array2 + sizeof(m_wd_array2)/sizeof(double));
static vector< vector<double> > table_Mwd_time;
static vector< vector<double> > table_Mwd_temp;
static vector< vector<double> > table_Mwd_u;
static vector< vector<double> > table_Mwd_g;
static vector< vector<double> > table_Mwd_r;
static vector< vector<double> > table_Mwd_i;
static vector< vector<double> > table_Mwd_z;
static vector< vector<double> > table_Mwd_Mbol;

void read_global_table(string filename, vector< vector<double> > &table){

    ifstream ifs(filename.c_str());
    if (!ifs) cerr << "error: couldn't open file "<< filename <<endl;
    
    while (!ifs.eof()){
        string line;
        getline(ifs, line);
        //   cout<<line<<endl;
        istringstream iss(line);
        vector<double> row;
        copy(istream_iterator<double>(iss), istream_iterator<double>(), back_inserter(row));
        table.push_back(row);
    };    
    ifs.close();
}


//i_wd table of lower limit mass
//i_wd+1 table of upper limit mass
//j_wd lower limit mass, lower limit time
//j_wd+1 lower limit mass, upper limit time
//k_wd upper limit mass, lower limit time
//k_wd+1 upper limit mass, upper limit time
void find_nearest_neighbours_3d(int& ii_1, int& ii_2, int& jj_1, int& jj_2, int& kk_1, int& kk_2, 
                                    const double mass, const double time, 
                                    const vector<double>& m_array, const vector< vector<double> > &table_time){
     
//    cerr<< "find_nearest_neighbours_3d" <<endl;                                   
    int i_1, j_1, k_1; 
    int i_2, j_2, k_2; 
    if (ii_1 < 0){
        // first iteration  
        
        int length_mass = m_array.size();
        if (mass <= m_array[0]){
             i_1 = 0;       
             i_2 = 0;       
        }
        else if (mass >= m_array[length_mass-1]){
             i_1 = length_mass-1;
             i_2 = length_mass-1;             
        }
        else{        
            for (i_2 = 0; i_2 < length_mass; ++i_2){
                if (mass < m_array[i_2]) break;    
            }
            i_1 = i_2-1;
        }

        int length = table_time[i_1].size();
        if (time <= table_time[i_1][0]){
            j_1 = 0;
            j_2 = 0;
        }
        else if (time >= table_time[i_1][length-1]){
            j_1 = length-1;
            j_2 = length-1;
        }
        else{
            for (j_2 = 0; j_2 < length; ++j_2){
                if (time < table_time[i_1][j_2]) break;                 
            }        
            j_1 = j_2-1;
        }
        
        int length2 = table_time[i_2].size(); 
        if (time <= table_time[i_2][0]){
             k_1 = 0;
             k_2 = 0;
        }
        else if (time >= table_time[i_2][length2-1]){ 
            k_1 = length2 -1;
            k_2 = length2 -1;
        }
        else{
            for (k_2 = 0; k_2 < length2; ++k_2){
                if (time < table_time[i_2][k_2]) break;                 
            }     
            k_1 = k_2-1;
        } 
    }
    else{// after first iteration
        i_1 = ii_1;     
        i_2 = ii_2;     
        int length = table_time[i_1].size();        
        if (time <= table_time[i_1][0]){
             j_1 = 0;
             j_2 = 0;
        }
        else if(time >= table_time[i_1][length-1]){
             j_1 = length-1;
             j_2 = length-1;
        }
        else{
            for (j_2 = jj_1; j_2 < length; ++j_2){
                if (time < table_time[i_1][j_2]) break;                 
            }     
            j_1 = j_2-1;
        }
        
        int length2 = table_time[i_2].size();                
        if (time <= table_time[i_2][0]){
             k_1 = 0;
             k_2 = 0;             
        }
        else if(time >= table_time[i_2][length2-1]){
             k_1 = length2-1;
             k_2 = length2-1;
        }
        else{
            for (k_2 = kk_1; k_2 < length2; ++k_2){
                if (time < table_time[i_2][k_2]) break;                 
            }     
            k_1 = k_2-1;           
        }
    }
    
    ii_1 = i_1;
    ii_2 = i_2;
    jj_1 = j_1;
    jj_2 = j_2;
    kk_1 = k_1;
    kk_2 = k_2;
}

double interpol_irregular_grid(const int ii_1, const int ii_2, const int jj_1, 
                                const int jj_2, const int kk_1, const int kk_2, 
                                const double mass, const double time, 
                                const vector<double>& m_array, vector< vector<double> > &table_time, 
                                vector< vector<double> > &table_interpol){

//    cerr<< "interpol_irregular_grid" <<endl;                                   
    double grid_1, grid_2;  
    double grid_11 = table_interpol[ii_1][jj_1];
    double time_11 = table_time[ii_1][jj_1];  
    if (jj_1 != jj_2){     
        double grid_12 = table_interpol[ii_1][jj_2];  
        double time_12 = table_time[ii_1][jj_2];   
        grid_1 = grid_11 + (time - time_11)/(time_12 - time_11) * (grid_12 - grid_11);
    }
    else grid_1 = grid_11;
    
    if (ii_1 == ii_2){
        return grid_1;
    }

    double grid_21 = table_interpol[ii_2][kk_1];  
    double time_21 = table_time[ii_2][kk_1];   
    if (kk_1 != kk_2){
        double grid_22 = table_interpol[ii_2][kk_2];  
        double time_22 = table_time[ii_2][kk_2];  
        grid_2 = grid_21 + (time - time_21)/(time_22 - time_21) * (grid_22 - grid_21);   
    }
    else grid_2 = grid_21;
    
    double grid_interpol = grid_1 + (mass - m_array[ii_1])/(m_array[ii_2] - m_array[ii_1]) * (grid_2 - grid_1);   
    return grid_interpol; 
    
}              


//Rappaport, Verbunt, Joss 1983
// P**y = Pi**y + y*xi*(t-ti)
// t = ti + (Pmax**y - Pi**y)/y/xi
double P_next(double Mwd1, double Mwd2, double Pi, double ti, double t, double z) {

    double P = 0;
    Pi *= 8.64e4; // in sec
    double yr_to_sec = 3.1536E7;
    
    double xi_gw = -3.68e-6 * Mwd1 * Mwd2 / pow(Mwd1 + Mwd2, ONE_THIRD);
    double y_gw = 8./3.; 
    P = pow(pow(Pi,y_gw) + (t-ti)*yr_to_sec*y_gw*xi_gw, 1/y_gw);
    
    P /= 8.64e4;
    return P;

}

double start_time(double Mwd1, double Mwd2, double Pmax, double Pi, double ti, double z) {

    double t = 0;
    Pi *= 8.64e4;
    Pmax *= 8.64e4;
    double yr_to_sec = 3.1536E7;
    
    double xi_gw = -3.68e-6 * Mwd1 * Mwd2 / pow(Mwd1 + Mwd2, ONE_THIRD);
    double y_gw = 8./3.; 
    t = ti + (pow(Pmax,y_gw)-pow(Pi,y_gw)) /yr_to_sec/y_gw/xi_gw; 
    
    return t;
}


double period_RLOF(double Mwd1, double Mwd2, double Rwd1, double Rwd2) {

    double q1_3   = pow(Mwd2/Mwd1, ONE_THIRD);
    double q2_3   = q1_3*q1_3;  
    double ap_RLOF = Rwd1*(0.6*(1./q2_3) + log(1 + (1./q1_3)))/(0.49*(1./q2_3));
    double as_RLOF = Rwd2*(0.6*q2_3 + log(1 + q1_3))/(0.49*q2_3);
    double a_RLOF  = max(ap_RLOF,as_RLOF);
    double P = 0.116 * sqrt(pow(a_RLOF, 3.)/(Mwd1 + Mwd2)); // in days   
        
    return P;

}


int main(int argc, char* argv[]) {    
    //assuming file contains 250000 binaries, and a binary fraction of 50%, otherwise change coeff    
    //assuming a SeBa run with: ./SeBa -f 4 -m 0.95 -M 10 -> needed for birth_rate_coefficient in myutil.C 
    char filename[] = "/Users/silviato/Development/process_data/data_run_SeBa/data_r105/SeBa_r105_ag_wdwd.data";
//    double coeff= birth_rate_coefficient(false, 4, 0, 500000, 0.5);    
    
    // if M_wd<0.2, the colors of a 0.2Msun are taken, similarly to a WD mass more massive then 1.2
    double g_lim = 23; 
    double P_max = 100; // in days
    double z = SOLAR_METALICITY; 
    int set_resolution = 0; // desired resolution
        // 0: default
        // 1: 10*MW resolution, P_max = 100
        // 2: 50*MW resolution, all periods, distance < 200 pc

    for (int i = 1; i < argc; ++i) { 
        if (i + 1 != argc){
            if ( strcmp(argv[i], "-f") == 0) {
                strcpy(filename, argv[i + 1]);
                i++;
            } else if ( strcmp(argv[i], "-g") == 0) {
                g_lim = atof(argv[i + 1]);
                i++;
            } else if ( strcmp(argv[i], "-P") == 0) {
                P_max = atof(argv[i + 1]);
                i++;
            } else if (strcmp(argv[i], "-R") == 0) {
                set_resolution = atoi(argv[i + 1]);
                i++;                
            } else {
                cout << "Not enough or invalid arguments, please try again.\n";
                exit(0);
            }
        }
    }
    
    read_global_table("Bergeron_time.txt", table_Mwd_time);
    read_global_table("Bergeron_temp.txt", table_Mwd_temp);
    read_global_table("Bergeron_u.txt", table_Mwd_u);
    read_global_table("Bergeron_g.txt", table_Mwd_g);
    read_global_table("Bergeron_r.txt", table_Mwd_r);
    read_global_table("Bergeron_i.txt", table_Mwd_i);
    read_global_table("Bergeron_z.txt", table_Mwd_z);
    read_global_table("Bergeron_Mbol.txt", table_Mwd_Mbol);

    double coeff;
    if (set_resolution == 0){
        coeff= birth_rate_coefficient(false, 4, 0, 500000, 0.5);
    } else if (set_resolution == 1) {//high resolution
        coeff= birth_rate_coefficient(false, 4, 0, 0.1*500000, 0.5);
    } else if (set_resolution == 2) {//high resolution & d<0.2
        coeff= birth_rate_coefficient(false, 4, 0, 0.02*500000, 0.5);
        P_max = 1e10;
    } else {
        cout << "Resolution not defined, please try again.\n";
        exit(0);        
    }       

    
    cerr<<"# The g band magnitude limit is: "<<g_lim<<endl;
    cerr<<"# The maximum period is: "<<P_max<<"d"<<endl;
    cerr<<"# Data from: " << filename <<endl;
    cerr<<"# Resolution is: " << set_resolution << endl;
        
    double T_end = 1.35e10; //in yrs
    double semi =0, ecc = 0, formation_time_wd1 =0, formation_time_wd2 =0;
    double l, b, d = 0.;    

    double wd1_mass = 0, wd2_mass = 0;
    int wd1_id = 0, wd2_id = 0, bin_id = 0;
    int wd1_ii1 = 0, wd1_jj1 = 0, wd1_kk1 = 0;
    int wd1_ii2 = 0, wd1_jj2 = 0, wd1_kk2 = 0;
    int wd2_ii1 = 0, wd2_jj1 = 0, wd2_kk1 = 0;
    int wd2_ii2 = 0, wd2_jj2 = 0, wd2_kk2 = 0;

    int num = 0;   
    double dt;//to read in data
    
    ifstream is(filename, ios::in);
    if (!is) cerr << "error: couldn't open file "<< filename <<endl;
    
    do {
        is >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt; 
        is >> dt >> dt >> dt >> formation_time_wd1 >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt >> dt; 
        
        is >> bin_id >> dt >> dt >> formation_time_wd2 >> semi >> ecc >> 
                dt >> wd1_id >> wd1_mass >> dt >> dt >> dt >> 
                dt >> wd2_id >> wd2_mass >> dt >> dt >> dt; 
        
        wd1_ii1 = -1, wd1_jj1 = -1, wd1_kk1 = -1;
        wd1_ii2 = -1, wd1_jj2 = -1, wd1_kk2 = -1;
        wd2_ii1 = -1, wd2_jj1 = -1, wd2_kk1 = -1;
        wd2_ii2 = -1, wd2_jj2 = -1, wd2_kk2 = -1;

        if (wd1_id < 12 || wd1_id > 14 || wd2_id < 12 || wd2_id > 14) {
        
            cerr << "ERROR something wrong!!!" << endl;
            PRC(wd1_id);PRL(wd1_mass);
            PRC(wd2_id);PRL(wd2_mass);
            PRC(formation_time_wd1);PRC(formation_time_wd2);PRL(semi);
            exit(1);
        }
    
        double P = 0.116 * sqrt(pow(semi, 3.)/(wd1_mass + wd2_mass)); // in days 
        double Pi = P;
        double wd1_rad = wd_radius_PPE(wd1_mass);
        double wd2_rad = wd_radius_PPE(wd2_mass);

        double t_form = formation_time_wd2 * 1.e6; // formation of double white dwarf
        double t_form_wd1 = formation_time_wd1 * 1.e6; // formation of first white dwarf
        double t = t_form;
        if (Pi > P_max) t = start_time(wd1_mass, wd2_mass, P_max, Pi, t_form, z); 

        do {
            
            double next_random_time;
            double t_next = next_birth_interval_BP(t, 
            		     T_end, coeff, &next_random_time);
            num++;

            double t_gal = T_end - next_random_time;
            
            // position in the galaxy  
            random_position_BP(t_gal, &l, &b, &d);
            double A_V = get_A_V(l, b, d);
            double logd = log10(d);           

            if (set_resolution <2 || (set_resolution ==2 & d<0.2)){

            	if (next_random_time <= T_end) {
      
                    P = P_next(wd1_mass, wd2_mass, Pi, t_form, next_random_time, z); 
                    if (P < P_max) { 
    
                        // RLOF
                        double P_RLOF = period_RLOF(wd1_mass, wd2_mass, wd1_rad, wd2_rad);
                        if (P <= P_RLOF) break; 

//                    double t_gal = T_end - next_random_time;
//                    
//                    // position in the galaxy  
//                    random_position_BP(t_gal, &l, &b, &d);
//                    double A_V = get_A_V(l, b, d);
//                    double logd = log10(d);           
                                            
                        // WD colors  
                        double u_wd1, g_wd1, r_wd1, i_wd1, z_wd1, T_wd1, Mbol_wd1;
                        find_nearest_neighbours_3d(wd1_ii1, wd1_ii2, wd1_jj1, wd1_jj2, wd1_kk1, wd1_kk2, 
                                wd1_mass, next_random_time - t_form_wd1, m_wd_array, table_Mwd_time);
                        u_wd1 = interpol_irregular_grid(wd1_ii1, wd1_ii2, wd1_jj1, wd1_jj2, wd1_kk1, wd1_kk2, 
                                wd1_mass, next_random_time-t_form_wd1, m_wd_array, table_Mwd_time, table_Mwd_u);              
                        g_wd1 = interpol_irregular_grid(wd1_ii1, wd1_ii2, wd1_jj1, wd1_jj2, wd1_kk1, wd1_kk2, 
                                wd1_mass, next_random_time-t_form_wd1, m_wd_array, table_Mwd_time, table_Mwd_g);              
                        r_wd1 = interpol_irregular_grid(wd1_ii1, wd1_ii2, wd1_jj1, wd1_jj2, wd1_kk1, wd1_kk2, 
                                wd1_mass, next_random_time-t_form_wd1, m_wd_array, table_Mwd_time, table_Mwd_r);              
                        i_wd1 = interpol_irregular_grid(wd1_ii1, wd1_ii2, wd1_jj1, wd1_jj2, wd1_kk1, wd1_kk2, 
                                wd1_mass, next_random_time-t_form_wd1, m_wd_array, table_Mwd_time, table_Mwd_i);              
                        z_wd1 = interpol_irregular_grid(wd1_ii1, wd1_ii2, wd1_jj1, wd1_jj2, wd1_kk1, wd1_kk2, 
                                wd1_mass, next_random_time-t_form_wd1, m_wd_array, table_Mwd_time, table_Mwd_z);              
                        T_wd1 = interpol_irregular_grid(wd1_ii1, wd1_ii2, wd1_jj1, wd1_jj2, wd1_kk1, wd1_kk2, 
                                wd1_mass, next_random_time-t_form_wd1, m_wd_array, table_Mwd_time, table_Mwd_temp);             
                        Mbol_wd1 = interpol_irregular_grid(wd1_ii1, wd1_ii2, wd1_jj1, wd1_jj2, wd1_kk1, wd1_kk2, 
                               wd1_mass, next_random_time-t_form_wd1, m_wd_array, table_Mwd_time, table_Mwd_Mbol);             
 
                        //   d in kpc
                        u_wd1 += 10. + 5*logd + 1.58*A_V;
                        g_wd1 += 10. + 5*logd + 1.16*A_V;
                        r_wd1 += 10. + 5*logd + 0.84*A_V;
                        i_wd1 += 10. + 5*logd + 0.64*A_V;
                        z_wd1 += 10. + 5*logd + 0.45*A_V;
    
    
                        double u_wd2, g_wd2, r_wd2, i_wd2, z_wd2, T_wd2, Mbol_wd2;
                        find_nearest_neighbours_3d(wd2_ii1, wd2_ii2, wd2_jj1, wd2_jj2, wd2_kk1, wd2_kk2, 
                                wd2_mass, next_random_time - t_form, m_wd_array, table_Mwd_time);
                        u_wd2 = interpol_irregular_grid(wd2_ii1, wd2_ii2, wd2_jj1, wd2_jj2, wd2_kk1, wd2_kk2, 
                                wd2_mass, next_random_time-t_form, m_wd_array, table_Mwd_time, table_Mwd_u);              
                        g_wd2 = interpol_irregular_grid(wd2_ii1, wd2_ii2, wd2_jj1, wd2_jj2, wd2_kk1, wd2_kk2, 
                                wd2_mass, next_random_time-t_form, m_wd_array, table_Mwd_time, table_Mwd_g);              
                        r_wd2 = interpol_irregular_grid(wd2_ii1, wd2_ii2, wd2_jj1, wd2_jj2, wd2_kk1, wd2_kk2, 
                                wd2_mass, next_random_time-t_form, m_wd_array, table_Mwd_time, table_Mwd_r);              
                        i_wd2 = interpol_irregular_grid(wd2_ii1, wd2_ii2, wd2_jj1, wd2_jj2, wd2_kk1, wd2_kk2, 
                                wd2_mass, next_random_time-t_form, m_wd_array, table_Mwd_time, table_Mwd_i);              
                        z_wd2 = interpol_irregular_grid(wd2_ii1, wd2_ii2, wd2_jj1, wd2_jj2, wd2_kk1, wd2_kk2, 
                                wd2_mass, next_random_time-t_form, m_wd_array, table_Mwd_time, table_Mwd_z);              
                        T_wd2 = interpol_irregular_grid(wd2_ii1, wd2_ii2, wd2_jj1, wd2_jj2, wd2_kk1, wd2_kk2, 
                                wd2_mass, next_random_time-t_form, m_wd_array, table_Mwd_time, table_Mwd_temp);             
                        Mbol_wd2 = interpol_irregular_grid(wd2_ii1, wd2_ii2, wd2_jj1, wd2_jj2, wd2_kk1, wd2_kk2, 
                                wd2_mass, next_random_time-t_form, m_wd_array, table_Mwd_time, table_Mwd_Mbol);             
    
                        //   d in kpc
                        u_wd2 += 10. + 5*logd + 1.58*A_V;
                        g_wd2 += 10. + 5*logd + 1.16*A_V;
                        r_wd2 += 10. + 5*logd + 0.84*A_V;
                        i_wd2 += 10. + 5*logd + 0.64*A_V;
                        z_wd2 += 10. + 5*logd + 0.45*A_V;
    
                        
                        if (g_wd1 < g_lim || g_wd2 < g_lim) {                            
                            cout << bin_id << " " << P << " " 
                            << l << " " << b << " " << d << " " 
                            << A_V << " "
                            << wd1_mass << " "  << wd1_rad << " " << T_wd1 << " "  
                            << wd2_mass << " "  << wd2_rad << " " << T_wd2 << " " 
                            << u_wd1 << " " << g_wd1 << " "  << r_wd1 << " " 
                            << i_wd1 << " " << z_wd1 << " " 
                            << u_wd2  << " " << g_wd2  << " " << r_wd2  << " "
                            << i_wd2  << " " << z_wd2  << " "
                            << Mbol_wd1 << " " << Mbol_wd2 << " " 
                            << endl;
                        } // cout
                    } // P < P_max
                } // next_random_time <= T_end
            } // resolution = 0,1 or (resolution = 2 & d<200 pc)  
            // evolve to end of "next birth interval"
            t = t_next;
            P = P_next(wd1_mass, wd2_mass, Pi, t_form, t, z);
        } while (t < T_end); // do loop

    } while (!is.eof()); // end of file
    is.close();    
    cout << "# Total = " << num << endl;
}
