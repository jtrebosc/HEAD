#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <vector>
#include <gsl/gsl_multimin.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

#include "simpson_rw.hpp"

// 1H Larmor Base Frequency
const double HLarmor = 1.;

struct DataStruct{
    int TD1;
    int TD2;
    int hetTD1;
    int hetTD2;
    double sw;
    double left;
    double right;
    double sw1;
    double left1;
    double right1;
    double hetsw1;
    double hetleft;
    double hetright;
    double hetleft1;
    double hetright1;
    double lambda;
    vector<int> index;
    vector< vector<double> > spectrum;
    vector< vector<double> > spectrum_scaled;
    vector<double> F2_sum;
    vector<double> F1_sum;

    vector< vector<double> > hetspectrum;
    vector<double> hetF1_sum;
    vector<int> hetindex;

    int slice;
} ;

DataStruct importSimpson(char* filename, char *filename1, int datatype) {

    DataStruct spec;
    SSpectrum sspec;
    std:string sfilename;

    // read the HEAD spectrum
    sfilename = filename;
    if (! readSimpsonFile(sfilename, sspec)) {
        cout << "Cannot read input file: " << sfilename;
        exit(1);
    }

    // initialize spec structure elements
    spec.lambda = 0.;
    spec.TD2 = sspec.NP;
    spec.sw = sspec.SW;
    spec.left = -sspec.SW/2-sspec.REF;
    spec.right = sspec.SW/2*(sspec.NP-1)/(sspec.NP)-sspec.REF;
    spec.TD1 = sspec.NI;
    spec.left1 = -sspec.SW1/2-sspec.REF1;
    spec.right1 = sspec.SW1/2*(sspec.NI-1)/(sspec.NI)-sspec.REF1;
    spec.sw1 = sspec.SW1;
    spec.F2_sum.resize(spec.TD2,0.);
    spec.F1_sum.resize(spec.TD2,0.);
    spec.spectrum.resize(spec.TD2);
    for(int i=0;i<spec.TD2;i++){
        spec.spectrum[i].resize(spec.TD2,0.);
    }
    spec.spectrum_scaled.resize(spec.TD2);
    for(int i=0;i<spec.TD2;i++){
        spec.spectrum_scaled[i].resize(spec.TD2,0.);
    }

    // store spec.spectrum from matrix

    for(int j=0;j<spec.TD1;j++){
        for(int i=0; i<spec.TD2; i++){
                int F1_index=-spec.TD1/2+i+j;
                if((F1_index>0)&&(F1_index<spec.TD2)){
                    spec.spectrum[F1_index][i] = sspec.ComplexData[j][i].real();
                    if(spec.spectrum[F1_index][i]<0.)
                        spec.spectrum[F1_index][i]=0.;
                }
        }
    }

    // -------------  correlation spectrum  -----------
    // DataStruct elements to update
    // vector< vector<double> > hetspectrum;
    // vector<double> hetF1_sum;
    // vector<int> hetindex;
    // double hetsw1;
    // double hetleft1;
    // double hetright1;

    // read the corretion spectrum
    // SSpectrum sspec1;
    sspec.ComplexData.clear();
    sfilename = filename1;
    if (! readSimpsonFile(sfilename, sspec)) {
        cout << "Cannot read input file: " << sfilename;
        exit(1);
    }
    spec.hetTD2 = sspec.NP;
    spec.hetTD1 = sspec.NI;
    spec.hetleft1 = -sspec.SW1/2-sspec.REF1;
    spec.hetright1 = sspec.SW1/2*(sspec.NI-1)/(sspec.NI)-sspec.REF1;
    spec.hetsw1 = sspec.SW1;

    if ((spec.TD2 != spec.hetTD2) && (spec.left != spec.hetleft) && (spec.right != spec.hetright)) {
        FILE *error_file;
        error_file=fopen("error.txt","a");
        fprintf(error_file, "\nERROR: the correlation spectrum must have same F2 acquisition window as HEAD spectrum.\n");
        fclose(error_file);
        exit(1);
    }


    if (datatype==2) {
        if ((spec.hetTD1 != spec.hetTD2) && (spec.hetleft1 != spec.hetleft) && (spec.hetright1 != spec.hetright)) {
            FILE *error_file;
            error_file=fopen("error.txt","a");
            fprintf(error_file, "\nERROR: the HOMCOR must have same acquisition window in F1 and F2.\n");
            fclose(error_file);
            exit(1);
        }
    }
    if (datatype==3) {
        if (spec.hetTD1 != spec.hetTD2*2) {
            FILE *error_file;
            error_file=fopen("error.txt","a");
            fprintf(error_file, "\nERROR: for DQSQ, the F1 dimension must have twice the point of F2 dimension.\n");
            fclose(error_file);
            exit(1);
        }
    }

    int TD1X = spec.hetTD1;
    if(datatype==3)//DSSQ into a SQSQ format
        TD1X/=2;
    spec.hetspectrum.resize(TD1X);
    for(int i=0;i<TD1X;i++){
        spec.hetspectrum[i].resize(spec.TD2,0.);
    }
    spec.hetF1_sum.resize(TD1X,0.);

    // store spec.hetspectrum from matrix1
    if (datatype!=3) {
        for(int j=0;j<TD1X;j++){
            for(int i=0; i<spec.TD2; i++){
                spec.hetspectrum[j][i] = sspec.ComplexData[j][i].real();
            }
        }
    } else { //DQSQ reading into SQSQ
        for(int j=0;j<TD1X*2;j++){
            for(int i=0; i<spec.TD2; i++){
                int F1_index=-i+j;
                if((F1_index>0)&&(F1_index<spec.TD2)){
                    spec.hetspectrum[F1_index][i] = sspec.ComplexData[j][i].real();
                    if(spec.hetspectrum[F1_index][i]<0.)
                        spec.hetspectrum[F1_index][i]=0.;
                }
            }
        }
    }
    return spec ;
}


vector<vector<double>> read_csv_matrix(char *filename, vector<double> *x_values, vector<double> *y_values) {
    // Reads a csv spectrum generated with ssNake->export->CSV
    // The csv data is transposed meaning rows are D2 and columns D1

    // Open the file
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open the file '" << filename << "'." << endl;
        exit(1);
    }

    vector<vector<double>> matrix;

    string line;

    // Read the first line (y-axis values)
    if (getline(file, line)) {
        istringstream iss(line);
        string token;
        if (getline(iss, token, ',')) {
            if (token.find("#") == 0) {
                y_values->push_back(stod(token.erase(0, 1)));
                while (getline(iss, token, ',')) {
                    y_values->push_back(stod(token));
                }
            }
        }
    }

    // Read the rest of the file skipping imaginary part
    while (getline(file, line)) {
        istringstream iss(line);
        string token;
        vector<double> row;

        // Read x-axis value (stored as first column value)
        if (getline(iss, token, ',')) {
            x_values->push_back(stod(token));
        }

        // Read complex values in the row
        int i = 0;
        while (getline(iss, token, ',')) {
            double real = stod(token);
            if (getline(iss, token, ',')) {
                double imag = stod(token);
            }
            row.emplace_back(real);
        }

        matrix.push_back(row);
    }
    file.close();
    return matrix;
}

DataStruct readMatrixFile2 (char* filename, char* filename1, int datatype) {
    // readMatrixFile2 reads file  stored in CSV format as exported by ssNake:
    // X values are stored as 1st column
    // Y values are stored as 1st row starting with #
    // Arguments:
    // filename: name of HEAD file (echo decay 2D spctrum)
    // left: value of left most point along X axis of filename
    // right: value of right most point along X axis of filename
    // filename1: name of HETCOR, HOMCOR or DQSQ file whose resolution will be improved
    // left1: value of left most point along Y axis of filename1 (het/hom/DQ)correlation spectrum
    // right1: value of right most point along Y axis of filename1 (het/hom/DQ)correlation spectrum
    // filename1 X axis should be the same as HEAD spectrum (echo decay)
    vector<double> x_values, y_values;
    vector<double> x1_values, y1_values;
    vector<vector<double>> matrix;
    vector<vector<double>> matrix1;

    matrix = read_csv_matrix(filename, &x_values, &y_values);
    double left = x_values[0];
    double right = x_values[x_values.size()-1];
    double left1 = y_values[0];
    double right1 = y_values[y_values.size()-1];

    matrix1 = read_csv_matrix(filename1, &x1_values, &y1_values);
    double hetleft = x1_values[0];
    double hetright = x1_values[x_values.size()-1];
    double hetleft1 = y1_values[0];
    double hetright1 = y1_values[y1_values.size()-1];

    int TD2 = x_values.size();
    int TD1 = y_values.size();
    int TD1X = y1_values.size();
    if(datatype==3)//DSSQ into a SQSQ format
        TD1X/=2;

    // Create DataStruct and update sizes
    DataStruct spec;
    spec.lambda = 0.;
    spec.TD2=TD2;
    spec.TD1=TD1;
    spec.hetTD1=TD1X;
    spec.left = right<left?right:left;
    spec.right = right>left?right:left;
    spec.left1 = right1<left1?right1:left1;
    spec.right1 = right1>left1?right1:left1;
    spec.hetleft = hetright<hetleft?hetright:hetleft;
    spec.hetright = hetright>hetleft?hetright:hetleft;
    spec.hetleft1 = hetright1<hetleft1?hetright1:hetleft1;
    spec.hetright1 = hetright1>hetleft1?hetright1:hetleft1;
    spec.sw = (right-left)/(TD2-1)*TD2;
    spec.sw1 = (right1-left1)/(TD1-1)*TD1;
    spec.hetsw1 = (hetright1-hetleft1)/(TD1X-1)*TD1X;
    spec.F2_sum.resize(TD2,0.);
    spec.F1_sum.resize(TD2,0.);

    if ((spec.TD2 != spec.hetTD2) && (spec.left != spec.hetleft) && (spec.right != spec.hetright)) {
        FILE *error_file;
        error_file=fopen("error.txt","a");
        fprintf(error_file, "\nERROR: the correlation spectrum must have same F2 acquisition window as HEAD spectrum.\n");
        fclose(error_file);
        exit(1);
    }


    if (datatype==2) {
        if ((spec.hetTD1 != spec.hetTD2) && (spec.hetleft1 != spec.hetleft) && (spec.hetright1 != spec.hetright)) {
            FILE *error_file;
            error_file=fopen("error.txt","a");
            fprintf(error_file, "\nERROR: the HOMCOR must have same acquisition window in F1 and F2.\n");
            fclose(error_file);
            exit(1);
        }
    }
    if (datatype==3) {
        if (spec.hetTD1 != spec.hetTD2*2) {
            FILE *error_file;
            error_file=fopen("error.txt","a");
            fprintf(error_file, "\nERROR: for DQSQ, the F1 dimension must have twice the point of F2 dimension.\n");
            fclose(error_file);
            exit(1);
        }
    }

    spec.spectrum.resize(TD2);
    spec.spectrum_scaled.resize(TD2);
    for(int i=0;i<TD2;i++){
        spec.spectrum[i].resize(TD2,0.);
        spec.spectrum_scaled[i].resize(TD2,0.);
    }
    //HETCOR spectrum to be fitted
    spec.hetspectrum.resize(TD1X);
    for(int i=0;i<TD1X;i++){
        spec.hetspectrum[i].resize(TD2,0.);
    }
    spec.hetF1_sum.resize(TD1X,0.);

    // store spec.spectrum from matrix
    for(int j=0;j<TD1;j++){
        for(int i=0; i<TD2; i++){
            int F1_index=-TD1/2+i+j;
            if((F1_index>0)&&(F1_index<TD2)){
                spec.spectrum[F1_index][i] = matrix[i][j];
                if(spec.spectrum[F1_index][i]<0.)
                    spec.spectrum[F1_index][i]=0.;
                spec.spectrum_scaled[F1_index][i]=spec.spectrum[F1_index][i];
            }
        }
    }

    // update spectrum_scaled
    for(int i=0;i<TD2;i++){
        double maximum=0;
        for(int j=0;j<TD2;j++){
            if(spec.spectrum_scaled[j][i]>maximum)
                maximum=spec.spectrum_scaled[j][i];
        }
        if(maximum>0.){
            for(int j=0;j<TD2;j++){
            spec.spectrum_scaled[j][i]= spec.spectrum_scaled[j][i]/maximum;
            }
        }
    }

    // store spec.hetspectrum from matrix1
    if (datatype!=3) {
        for(int j=0;j<TD1X;j++){
            for(int i=0; i<TD2; i++){
                spec.hetspectrum[j][i] = matrix1[i][j];
            }
        }
    } else { //DQSQ reading into SQSQ
        for(int j=0;j<TD1X*2;j++){
            for(int i=0; i<TD2; i++){
                int F1_index=-i+j;
                if((F1_index>0)&&(F1_index<TD2)){
                    spec.hetspectrum[F1_index][i] = matrix1[i][j];
                    if(spec.hetspectrum[F1_index][i]<0.)
                        spec.hetspectrum[F1_index][i]=0.;
                }
            }
        }
    }
    return spec ;

}

DataStruct read_totxt_file_meta(char *totxt_filename, char *totxt_filename1, int datatype) {
    FILE *fp;
    char buffer[256], word[24], pound;
    int i=0, j=0, TD1, TD2, TD1X, k, junk;
    double left, right, left1, right1, hetleft1, hetright1;

    //reading the header file to get TD1, TD2, and the spectral width and offset in F2.
    fp=fopen(totxt_filename,"r");
    if(fp==NULL){
        FILE *error_file;
        error_file=fopen("error.txt","a");
        fprintf(error_file, "\nERROR: totxt spectrum file '%s' not found\n", totxt_filename);
        fclose(error_file);
        exit(1);
    }

    int state=0;
    while ((fgets(buffer, sizeof(buffer), fp) != NULL)||(state<3)) {
        if(buffer[0]=='#'){
                sscanf(buffer,"%c %s",&pound,word);
                if(strcmp(word,"F2LEFT")==0){
                    sscanf(buffer,"%c %s %s %lf %s %s %s %lf",&pound,word,word,&left,word,word,word,&right);
                    sprintf(word,"void");
                    state++;
                }
                else if(strcmp(word,"F1LEFT")==0){
                    sscanf(buffer,"%c %s %s %lf %s %s %s %lf",&pound,word,word,&left1,word,word,word,&right1);
                    sprintf(word,"void");
                    state++;
                }
                else if(strcmp(word,"NROWS")==0){
                    sscanf(buffer,"%c %s %s %d",&pound,word,word,&TD1);
                    sprintf(word,"void");
                    state++;
                }
                else if(strcmp(word,"NCOLS")==0){
                    sscanf(buffer,"%c %s %s %d",&pound,word,word,&TD2);
                    sprintf(word,"void");
                    state++;
                }
            }
    }
    fclose(fp);

    // read second file header
    state=0;
    fp=fopen(totxt_filename1,"r");
    if(fp==NULL){
        FILE *error_file;
        error_file=fopen("error.txt","a");
        fprintf(error_file, "\nERROR: totxt spectrum file '%s' not found\n", totxt_filename1);
        fclose(error_file);
        exit(1);
    }
    while ((fgets(buffer, sizeof(buffer), fp) != NULL)||(state<3)) {
        if(buffer[0]=='#'){
                sscanf(buffer,"%c %s",&pound,word);
                if(strcmp(word,"F1LEFT")==0){
                    sscanf(buffer,"%c %s %s %lf %s %s %s %lf",&pound,word,word,&hetleft1,word,word,word,&hetright1);
                    sprintf(word,"void");
                    state++;
                }
                else if(strcmp(word,"NROWS")==0){
                    sscanf(buffer,"%c %s %s %d",&pound,word,word,&TD1X);
                    sprintf(word,"void");
                    state++;
                    if((datatype==2)&&(TD1X!=TD2)){
                        FILE *error_file;
                        error_file=fopen("error.txt","a");
                        fprintf(error_file, "\nERROR: the two TD2 values do not match\n");
                        fclose(error_file);
                        exit(1);
                    }
                    else if((datatype==3)&&(TD1X!=2*TD2)){
                        FILE *error_file;
                        error_file=fopen("error.txt","a");
                        fprintf(error_file, "\nERROR: the two TD2 values do not match\n");
                        fclose(error_file);
                        exit(1);
                    }
                }

                else if(strcmp(word,"NCOLS")==0){
                    sscanf(buffer,"%c %s %s %d",&pound,word,word,&junk);
                    sprintf(word,"void");
                    state++;
                    if(junk!=TD2){
                        FILE *error_file;
                        error_file=fopen("error.txt","a");
                        fprintf(error_file, "\nERROR: the two TD2 values do not match\n");
                        fclose(error_file);
                        exit(1);
                    }
                }
            }
    }
    fclose(fp);

    if(datatype==3)//DSSQ into a SQSQ format
        TD1X/=2;

    //creating the data structures
    DataStruct spec;
    spec.TD2=TD2;
    spec.TD1=TD1;
    spec.hetTD1=TD1X;
    spec.left = (right>left)?left:right;
    spec.right = (right<left)?left:right;
    spec.left1 = (right1>left1)?left1:right1;
    spec.right1 = (right1<left1)?left1:right1;
    spec.hetleft1 = (hetright1>hetleft1)?hetleft1:hetright1;
    spec.hetright1 = (hetright1<hetleft1)?hetleft1:hetright1;
    spec.spectrum.resize(TD2);
    spec.spectrum_scaled.resize(TD2);
    spec.F2_sum.resize(TD2,0.);
    spec.F1_sum.resize(TD2,0.);
    for(i=0;i<TD2;i++){
        spec.spectrum[i].resize(TD2,0.);
        spec.spectrum_scaled[i].resize(TD2,0.);
    }
    //HETCOR spectrum to be fitted
    spec.hetspectrum.resize(TD1X);
    for(i=0;i<TD1X;i++){
        spec.hetspectrum[i].resize(TD2,0.);
    }
    spec.hetF1_sum.resize(TD1X,0.);

    //reading the spectrum intensities and storing them as a sheared spectrum TD2xTD2
    fp=fopen(totxt_filename,"r");
    for(j=0;j<TD1;j++){
        for(i=0; i<TD2; i++){
            fgets(buffer, sizeof(buffer), fp);
            if(buffer[0]=='#'){
                i--;
            }
            else{
                int F1_index=-TD1/2+i+j;
                if((F1_index>0)&&(F1_index<TD2)){
                    sscanf(buffer,"%lf",&spec.spectrum[F1_index][i]);
                    if(spec.spectrum[F1_index][i]<0.)
                        spec.spectrum[F1_index][i]=0.;

                    spec.spectrum_scaled[F1_index][i]=spec.spectrum[F1_index][i];
                }
            }
        }
    }
    fclose(fp);
    for(i=0;i<TD2;i++){
        double maximum=0;
        for(j=0;j<TD2;j++){
            if(spec.spectrum_scaled[j][i]>maximum)
                maximum=spec.spectrum_scaled[j][i];
        }
        if(maximum>0.){
        for(j=0;j<TD2;j++){
            spec.spectrum_scaled[j][i]= spec.spectrum_scaled[j][i]/maximum;
        }}
    }

    if (datatype!=3) {
        //reading the HETCOR spectrum intensities
        fp=fopen(totxt_filename1,"r");
        for(j=0;j<TD1X;j++){
            for(i=0; i<TD2; i++){
                fgets(buffer, sizeof(buffer), fp);
                if (buffer[0]=='#') {
                    i--;
                } else {
                        sscanf(buffer,"%lf",&spec.hetspectrum[j][i]);
                        if(spec.hetspectrum[j][i]<0.)
                            spec.hetspectrum[j][i]=0.;
                      /*  for(k=0;k<TD2;k++){
                            spec.spectrum_scaled[j][k][i]*=spec.hetspectrum[j][i];
                        }*/
                }
            }
        }
        fclose(fp);
    } else { //DQSQ reading into SQSQ
        fp=fopen(totxt_filename1,"r");
        for(j=0;j<TD1X*2;j++){
            for(i=0; i<TD2; i++){

                int F1_index=-i+j;

                fgets(buffer, sizeof(buffer), fp);
                if(buffer[0]=='#'){
                    i--;
                }
                else if((F1_index>=0)&&(F1_index<TD2)){
                        sscanf(buffer,"%lf",&spec.hetspectrum[F1_index][i]);
                        if(spec.hetspectrum[F1_index][i]<0.)
                            spec.hetspectrum[F1_index][i]=0.;
                       /* for(k=0;k<TD2;k++){
                            spec.spectrum_scaled[F1_index][k][i]*=spec.hetspectrum[F1_index][i];
                        }*/
                }
            }
        }
        fclose(fp);
    }

    return spec;
}

double RMSD(const gsl_vector* weights, void* params){
    //Cost function used to calculate the offset between the experimental F2 spectrum
    //and the predicted isotropic spectrum in F1. Also includes Tikhonov weighting.

    //gathering the relevant parameters from the data structure
    DataStruct *spec = (DataStruct*) params;
    int TD2=spec->TD2,i,j;
    int TD1=spec->TD1;
    vector<double> F1_sum(TD2,0.);
    vector<double> corr_weights(TD2,0.);
    double MSD=0., gradient=0., Euc_norm=0.;
     vector<double> slice(TD2,0.);
    for(i=0;i<TD2;i++){
        slice[i]=spec->hetspectrum[spec->hetindex[spec->slice]][i];
    }

    //Derivative of the weights, used for Tikhonov
    for(i=1;i<spec->index.size();i++){
        gradient+=abs(gsl_vector_get(weights,i-1)-gsl_vector_get(weights,i));
    }
    gradient*=spec->lambda;

    //To avoid there being a larger lambda for the tails of the peak shapes
    gradient*=spec->hetF1_sum[spec->hetindex[spec->slice]];

    //Calculation of the F1 spectrum
    //looping over the basis spectra
    for(i=0;i<spec->index.size();i++){
        int ii=spec->index[i];
        int start_j= (ii-TD1/2)*(ii>=(TD1/2));
        int end_j= (ii+TD1/2)*((ii+TD1/2)<TD2)+(TD2-1)*((ii+TD1/2)>=TD2);
       // Euc_norm+= pow(gsl_vector_get(weights,i),2.);

        //looping over the data points of each basis spectrum
        for(j=start_j;j<=end_j;j++){
                F1_sum[j]=F1_sum[j] + spec->hetspectrum[spec->hetindex[spec->slice]][ii]*spec->spectrum_scaled[j][ii]*abs(gsl_vector_get(weights,i));
        }
    }

    //Calculation of the mean squared deviation between the calculated F1 spectrum and the F2 spectrum
    for(i=0;i<TD2;i++){
        MSD=MSD + pow(F1_sum[i]-slice[i],2.);
    }

    return sqrt(MSD)+gradient;
    //return sqrt(MSD)+spec->lambda*Euc_norm;
}

void gradient(const gsl_vector *var, void *params, gsl_vector *df){
    //Gradient of the RMSD, used for GSL gradient optimizers

    DataStruct *spec = (DataStruct *) params;
    int var_size=spec->index.size(),i;
    double cost_0 = RMSD(var, params), val[2];

    for(i=0;i<var_size;i++){
        val[0]=gsl_vector_get(var,i);
        val[1]=val[0]+0.0001;
        gsl_vector_set(var,i,val[1]);
        gsl_vector_set(df, i, (RMSD(var, params)-cost_0)/(val[1]-val[0]));
        gsl_vector_set(var,i,val[0]);
    }
}

void gradient_fdf (const gsl_vector *var, void *params,double *f,gsl_vector *df){
    //Simultaneously calculated the gradient and the value of the RMSD
    //For the GSL gradient optimizers

    *f = RMSD(var, params);
    DataStruct *spec = (DataStruct *) params;
    int var_size=spec->index.size(),i;
    double cost_0 = RMSD(var, params), val[2];

    for(i=0;i<var_size;i++){
        val[0]=gsl_vector_get(var,i);
        val[1]=val[0]+0.0001;
        gsl_vector_set(var,i,val[1]);
        gsl_vector_set(df, i, (RMSD(var, params)-cost_0)/(val[1]-val[0]));
        gsl_vector_set(var,i,val[0]);
    }
}

void calc_HETCOR(const gsl_vector* weights, void* params, vector< vector<double> >& iso_spec, int F1_index){
    //Function used to integrate over the HEAD spectrum to produce the isotropic spectrum with appropriate intensities

    //Gathering parameters from the data structure
    DataStruct *spec = (DataStruct*) params;
    int TD2=spec->TD2,i,j,ii;
    int TD1=spec->TD1;
    vector<double> corr_weights(TD2,0.);
    vector<double> F2_sum(TD2,0.);

    //Looping over the basis spectra
    for(i=0;i<spec->index.size();i++){
        ii=spec->index[i];
        int start_j= (ii-TD1/2)*(ii>=(TD1/2));
        int end_j= (ii+TD1/2)*((ii+TD1/2)<TD2)+(TD2-1)*((ii+TD1/2)>=TD2);

        //looping over the data points of each basis spectrum
        for(j=start_j;j<=end_j;j++){
            F2_sum[ii]=F2_sum[ii] + spec->hetspectrum[spec->hetindex[spec->slice]][ii]*spec->spectrum_scaled[j][ii]*abs(gsl_vector_get(weights,i));
        }

        //saving the isotropic spectrum and weights
        iso_spec[F1_index][ii]=F2_sum[ii];
    }
}

int main(int argc, char* argv[]) {
    //parameter declaration
    int TD2, TD1, TD1X, i=0, j=0,k,junk, datatype;
    char totxt_filename[128], totxt_filename1[128], buffer[256], word[24], pound, *type;
    double left, right, width,lambda, delta;
    double left1, right1, width1, delta1;
    FILE *fp;

    //Interface to gather the totxt filename and the lambda parameter for regularization
    if (argc < 5) { // min 4 args: COR_type input_type HEAD_filename COR_filename
        printf("Select the type of 2D spectrum you wish to Hahn-echo assisted deconvolute:\n");
        printf("1 1H-detected HETCOR\n");
        printf("2 SQ-SQ correlation\n");
        printf("3 SQ-DQ correlation\n");
        scanf("%d",&datatype);

        printf("What is the type of input file ?\n\
          Either 'totxt' for the 2D CS-lineshape correlation spectrum converted using totxt?\n\
          Or 'csv' for ssNake csv export.\n\
          Or 'spe' for ssNake simpson export.");
        scanf("%s",type);
        printf("What is the filename for the 2D CS-lineshape correlation spectrum?\n");
        printf("Note that the digital resolution in F1 and F2 needs to be identical.\n");
        scanf("%s",totxt_filename);
    } else {
        sscanf(argv[1], "%lf",&datatype);
        type = argv[2];
        strcpy(totxt_filename, argv[3]);
        strcpy(totxt_filename1, argv[4]);
    }


    if (argc < 6) { // min 5 args: COR_type input_type HEAD_filename COR_filename lambda
        printf("Lambda value for regularization (0 for pure least squares)\n");
        scanf("%lf",&lambda);
    } else {
        sscanf(argv[4], "%lf",&lambda);
    }

// ----------------------------------------------------------------------------------------
    DataStruct spec ;
    if (strncmp(type, "totxt", 5) == 0) {
        spec = read_totxt_file_meta(totxt_filename, totxt_filename1, datatype);
    } else if (strncmp(type, "csv", 3) == 0) {
        spec = readMatrixFile2(totxt_filename, totxt_filename1, datatype);
    } else if (strncmp(type, "spe", 3) == 0) {
        spec = importSimpson(totxt_filename, totxt_filename1, datatype);
    } else {
        printf("type could not be determined. Exit!");
        exit(1);
    }

    spec.lambda=lambda;
    TD2 = spec.TD2;
    TD1 = spec.TD1;
    TD1X = spec.hetTD1;

    width=spec.sw;
    // cout << totxt_filename << " : " << spec.left << ", " << spec.right << "-> width=" << width << "\n";
    delta=width/TD2;

    width1=spec.sw1;
    // cout << totxt_filename << " : " << spec.left1 << ", " << spec.right1 << "-> width=" << width1 << "\n";
    // cout << totxt_filename1 << " : " << spec.hetleft1 << ", " << spec.hetright1 << "-> width=" << spec.hetsw1 << "\n";
    delta1=width1/TD2;

    vector<double> F1_sum(TD2,0.);
    // cout << "OK2------------------------------------------\n";

    //Saving the original 2D spectrum as a *.spe file
    fp=fopen("original_HEAD.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*HLarmor,left*HLarmor,HLarmor,TD1,width1,left1,HLarmor);
    for(i=0;i<TD1;i++){
        for(j=0;j<TD2;j++){
            int F1_index=-TD1/2+i+j;
            if((F1_index>0)&&(F1_index<TD2))
                fprintf(fp,"%lf   0.0\n",spec.spectrum[F1_index][j]);
            else fprintf(fp,"0.0   0.0\n");
        }
    }
    fprintf(fp,"%s","END");
    fclose(fp);

    //Saving the sheared 2D spectrum as a *.spe file
    fp=fopen("sheared_HEAD.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",\
                         TD2,width*HLarmor,left*HLarmor,HLarmor,TD2,width*HLarmor,left*HLarmor,HLarmor);
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",spec.spectrum[i][j]);
        }
    }
    fprintf(fp,"%s","END");
    fclose(fp);

    //Saving the original HETCOR spectrum as a *.spe file
    if(datatype!=3){
    fp=fopen("original_correlation.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*HLarmor,left*HLarmor,HLarmor,TD1X,width1*150.0,left1*150.,150.0);
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",spec.hetspectrum[i][j]);
        }
    }
    fprintf(fp,"%s","END");
    fclose(fp);
    }

    else{
    fp=fopen("original_correlation.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,10.0*HLarmor,10.0*HLarmor,HLarmor,TD2*2,20.0*HLarmor,20.0*HLarmor,HLarmor);
    for(i=0;i<TD2*2;i++){
        for(j=0;j<TD2;j++){
            int F1_index=-j+i;
            if((F1_index>=0)&&(F1_index<TD2))
                fprintf(fp,"%lf   0.0\n",spec.hetspectrum[F1_index][j]);
            else fprintf(fp,"0   0.0\n");
        }}
    fprintf(fp,"%s","END");
    fclose(fp);
    }

    //normalizing the scaled HEAD spectrum
    double max_F2=0.;
    fill(spec.F2_sum.begin(), spec.F2_sum.end(), 0.);
    fill(F1_sum.begin(), F1_sum.end(), 0.);
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.F2_sum[j]=spec.F2_sum[j]+spec.spectrum_scaled[i][j];
            F1_sum[i]=F1_sum[i]+spec.spectrum_scaled[i][j];
            if(spec.F2_sum[j]>max_F2)
                max_F2=spec.F2_sum[j];
        }
    }

    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.spectrum_scaled[i][j]=spec.spectrum_scaled[i][j]/max_F2;
        }
        spec.F2_sum[i]=spec.F2_sum[i]/max_F2;
        F1_sum[i]=F1_sum[i]/max_F2;
    }


        //normalizing the F2 spectrum, (sum of rows)
    fill(spec.F2_sum.begin(), spec.F2_sum.end(), 0.);
    fill(F1_sum.begin(), F1_sum.end(), 0.);
    max_F2=0.;
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.F2_sum[j]=spec.F2_sum[j]+spec.spectrum[i][j];
            F1_sum[i]=F1_sum[i]+spec.spectrum[i][j];
            if(spec.F2_sum[j]>max_F2)
                max_F2=spec.F2_sum[j];
        }
    }

    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.spectrum[i][j]=spec.spectrum[i][j]/max_F2;
        }
        spec.F2_sum[i]=spec.F2_sum[i]/max_F2;
        F1_sum[i]=F1_sum[i]/max_F2;
    }


    //Normalizing the HETCOR
    max_F2=0.;
    double max_F1=0.;
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            if(spec.hetspectrum[i][j]>max_F2)
                max_F2=spec.hetspectrum[i][j];
    }}
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            spec.hetspectrum[i][j]=spec.hetspectrum[i][j]/max_F2;
            spec.hetF1_sum[i]+=spec.hetspectrum[i][j];
    }
        if(spec.hetF1_sum[i]>max_F1)
            max_F1=spec.hetF1_sum[i];
    }
    for(i=0;i<TD1X;i++){
        spec.hetF1_sum[i]/=max_F1;
    }


    //deciding which datapoint intensities are to be optimized
    //Here a default minimum intensity of 1% is used, this can be changed.
    for(i=0;i<TD2;i++){
         if(spec.F2_sum[i]>0.01){
            spec.index.push_back(i);
         }
    }
    for(i=0;i<TD1X;i++){
         if(spec.hetF1_sum[i]>0.05){
            spec.hetindex.push_back(i);
         }
    }

    vector< vector<double> > iso_spec;
    iso_spec.resize(TD1X);
    for(i=0;i<TD1X;i++){
        iso_spec[i].resize(TD2,0.);
    }

    int done=0;
    #pragma omp parallel for
    for(int J=0;J<spec.hetindex.size();J++){
    int I,K;
    DataStruct spec_temp=spec;
    printf("%d weights to optimize\n",spec.index.size());
    //setting the initial weights for the Tikhonov minimization to 1
    gsl_vector *weights;
    weights = gsl_vector_alloc (spec.index.size());
    spec_temp.slice=J;
    for(I=0;I<spec.index.size();I++){
        float wt=1.0;
        gsl_vector_set(weights,I,wt);
    }

    //The GSL gradient minimization is done here. This closely mirrors the example from
    //the GSL manual. BFGS, CG, and steepest descent algorithms all work well, but simplex
    //methods do not. If you chose to change the algorithm, you can uncomment the appropriate line below
    //The gsl_multimin_fdfminimizer_set parameters may need to be tweaked in that case.

    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
   // const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_steepest_descent;
    gsl_multimin_function_fdf minex_func; //declaring the minimizer
    minex_func.n = spec.index.size();  //number of variables
    minex_func.f = RMSD; //cost function assignment
    minex_func.df = gradient; //gradient calculation function
    minex_func.fdf = gradient_fdf; //simultaneous gradient and cost calculation function
    minex_func.params = &spec_temp; //parameter assignment (non-optimized things)
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, minex_func.n); //pointer for the minimizer
    gsl_multimin_fdfminimizer_set (s, &minex_func, weights, 0.1, 1e-0);  //Assigning the initial weights, functions, etc.
    double cost,cost_min=100000000.; //initial cost is set arbitrarily high

    size_t iter = 0;
    int status;

    //Optimization loop with a maximum of 1 million steps
    do{
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status)
            break;

        //maximum allowable gradient convergence criterion
        status = gsl_multimin_test_gradient (s->gradient, 0.00000001);

        if (status == GSL_SUCCESS){
            break;
        }
        cost=s->f;

        //Print to the screen if there is a better solution
          if(cost<cost_min){
            printf ("%d/%d %d %f \n",done,spec.hetindex.size(),iter,s->f);
            cost_min=cost;
        }
    }  while (status == GSL_CONTINUE && iter < 10000);

    printf("Converged!\n");
    done++;
    calc_HETCOR(s->x,&spec_temp,iso_spec,spec.hetindex[J]);
    gsl_vector_free(weights);
    gsl_multimin_fdfminimizer_free(s);
    }

    //Saving the isotropic spectrum as a *.spe file
    if(datatype!=3){
    fp=fopen("isotropic_correlation_F2only.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*HLarmor,left*HLarmor,HLarmor,TD1X,width1*HLarmor,left1*HLarmor,HLarmor);
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",iso_spec[i][j]);
        }
    }
    fprintf(fp,"%s","END");
    fclose(fp);
    }

    else{
    fp=fopen("isotropic_correlation_F2only.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,10.0*HLarmor,10.0*HLarmor,HLarmor,TD2*2,20.0*HLarmor,20.0*HLarmor,HLarmor);
    for(i=0;i<TD2*2;i++){
        for(j=0;j<TD2;j++){
            int F1_index=-j+i;
            if((F1_index>=0)&&(F1_index<TD2))
                fprintf(fp,"%lf   0.0\n",iso_spec[F1_index][j]);
            else fprintf(fp,"0   0.0\n");
        }}
    fprintf(fp,"%s","END");
    fclose(fp);
    }

    if(datatype==1)//HETCOR, onle F2 deconvolution
        return 0;

    //transpose
    for(i=0;i<TD2;i++){
        for(j=0;j<TD2;j++){
            spec.hetspectrum[i][j]=iso_spec[j][i];
        }
    }

    done=0;
    #pragma omp parallel for
    for(int J=0;J<spec.hetindex.size();J++){
    int I,K;
    DataStruct spec_temp=spec;
    printf("%d weights to optimize\n",spec.index.size());
    //setting the initial weights for the Tikhonov minimization to 1
    gsl_vector *weights;
    weights = gsl_vector_alloc (spec.index.size());
    spec_temp.slice=J;
    for(I=0;I<spec.index.size();I++){
        float wt=1.0;
        gsl_vector_set(weights,I,wt);
    }

    //The GSL gradient minimization is done here. This closely mirrors the example from
    //the GSL manual. BFGS, CG, and steepest descent algorithms all work well, but simplex
    //methods do not. If you chose to change the algorithm, you can uncomment the appropriate line below
    //The gsl_multimin_fdfminimizer_set parameters may need to be tweaked in that case.

    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
    //const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
   // const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_steepest_descent;
    gsl_multimin_function_fdf minex_func; //declaring the minimizer
    minex_func.n = spec.index.size();  //number of variables
    minex_func.f = RMSD; //cost function assignment
    minex_func.df = gradient; //gradient calculation function
    minex_func.fdf = gradient_fdf; //simultaneous gradient and cost calculation function
    minex_func.params = &spec_temp; //parameter assignment (non-optimized things)
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, minex_func.n); //pointer for the minimizer
    gsl_multimin_fdfminimizer_set (s, &minex_func, weights, 0.1, 1e-0);  //Assigning the initial weights, functions, etc.
    double cost,cost_min=100000000.; //initial cost is set arbitrarily high

    size_t iter = 0;
    int status;

    //Optimization loop with a maximum of 1 million steps
    do{
        iter++;
        status = gsl_multimin_fdfminimizer_iterate (s);

        if (status)
            break;

        //maximum allowable gradient convergence criterion
        status = gsl_multimin_test_gradient (s->gradient, 0.00000001);

        if (status == GSL_SUCCESS){
            break;
        }
        cost=s->f;

        //Print to the screen if there is a better solution
          if(cost<cost_min){
            printf ("%d/%d %d %f \n",done,spec.hetindex.size(),iter,s->f);
            cost_min=cost;
        }
    }  while (status == GSL_CONTINUE && iter < 10000);

    printf("Converged!\n");
    done++;
    calc_HETCOR(s->x,&spec_temp,iso_spec,spec.hetindex[J]);
    gsl_vector_free(weights);
    gsl_multimin_fdfminimizer_free(s);
    }

    //Saving the isotropic spectrum as a *.spe file
    if(datatype!=3){
    fp=fopen("isotropic_HOMCOR_spectrum.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,width*HLarmor,left*HLarmor,HLarmor,TD2,width*HLarmor,left*HLarmor,HLarmor);
    for(i=0;i<TD1X;i++){
        for(j=0;j<TD2;j++){
            fprintf(fp,"%lf   0.0\n",iso_spec[i][j]);
        }
    }
    fprintf(fp,"%s","END");
    fclose(fp);
    }

    else{
    fp=fopen("isotropic_HOMCOR_spectrum.spe","w");
    fprintf(fp,"SIMP\nNP=%d\nSW=%.2lf\nX0=%.2lf\nSF=%.2lf\nNI=%d\nSW1=%.2lf\nX0_F1=%.2lf\nSF1=%.2lf\nTYPE=SPE\nDATA\n",TD2,10.0*HLarmor,10.0*HLarmor,HLarmor,TD2*2,20.0*HLarmor,20.0*HLarmor,HLarmor);
    for(i=0;i<TD2*2;i++){
        for(j=0;j<TD2;j++){
            int F1_index=-j+i;
            if((F1_index>=0)&&(F1_index<TD2))
                fprintf(fp,"%lf   0.0\n",iso_spec[F1_index][j]);
            else fprintf(fp,"0   0.0\n");
        }}
    fprintf(fp,"%s","END");
    fclose(fp);
    }

    return 0;
}
