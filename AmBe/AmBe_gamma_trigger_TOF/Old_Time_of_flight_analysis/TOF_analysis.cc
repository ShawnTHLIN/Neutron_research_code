// Programmer: Dr. J (aka Jason Nattress)
// Purpose:
// Date: 2/3/2018
// Compliation command: 
// g++ mainBinary.cc WaveClass.cc HistoClass.cc `root-config --libs --cflags` -o main -O3 -ftree-vectorize
// Execution command: ./main /path/to/file outputFileName
// output file will be store in directory called "processed"
// to do: need to fix analyze class to mesh with new HistoClass.... done... I think
/////////////////////////////////////////////////////////
//**#include "analyze.h"
//**#include "HistoClass.h"
//#include "TrigEJ309.h"
//#include "VerifyEJ309.h"
//#include "MainEJ290.h"
//#include "constants.h"
//#include "main.h"
//#include <TH2.h>
//#include <TStyle.h>
//#include <TCanvas.h>
//#include <TCutG.h>

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <time.h>
#include <vector>
#include <thread>
//#include <mpi.h>

#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include "TCanvas.h"


using namespace std;

//24 hr 1.5 m 
//const string MainFileName =     "20211222_AmBe_400ns/UNFILTERED/Data_CH0@DT5730B_1173_20211221AmBe_TH_80_HV1550_D150_400ns_2.bin";
//const string ATriggerFileName = "20211222_AmBe_400ns/UNFILTERED/Data_CH4@DT5730B_1173_20211221AmBe_TH_80_HV1550_D150_400ns_2.bin";
//const string BTriggerFileName = "20211222_AmBe_400ns/UNFILTERED/Data_CH6@DT5730B_1173_20211221AmBe_TH_80_HV1550_D150_400ns_2.bin";

// 4 hr 1.5 m without coincidence
//const string MainFileName =     "20211222_AmBe_no_coin/UNFILTERED/Data_CH0@DT5730B_1173_20211221AmBe_TH_80_HV1550_D150_No_coin.bin";
//const string ATriggerFileName = "20211222_AmBe_no_coin/UNFILTERED/Data_CH4@DT5730B_1173_20211221AmBe_TH_80_HV1550_D150_No_coin.bin";
//const string BTriggerFileName = "20211222_AmBe_no_coin/UNFILTERED/Data_CH6@DT5730B_1173_20211221AmBe_TH_80_HV1550_D150_No_coin.bin";

// 2021/12/24  40 hrs D=1m d =3 cm
//const string MainFileName =     "20211224_AmBe_TH80_HV1550_d3cmD1m_400ns_40hr/UNFILTERED/Data_CH0.bin";
//const string ATriggerFileName = "20211224_AmBe_TH80_HV1550_d3cmD1m_400ns_40hr/UNFILTERED/Data_CH4.bin";
//const string BTriggerFileName = "20211224_AmBe_TH80_HV1550_d3cmD1m_400ns_40hr/UNFILTERED/Data_CH6.bin";

// 2021/12/26  9 hrs D=1m d =3 cm
//const string MainFileName =     "20211226_AmBe_TH80_HV1550_d3cmD1m_400ns_9hr/UNFILTERED/Data_CH0.bin";
//const string ATriggerFileName = "20211226_AmBe_TH80_HV1550_d3cmD1m_400ns_9hr/UNFILTERED/Data_CH4.bin";
//const string BTriggerFileName = "20211226_AmBe_TH80_HV1550_d3cmD1m_400ns_9hr/UNFILTERED/Data_CH6.bin";

// 2021/12/26  24 hrs D=1m d =5 cm
const string MainFileName =     "20211226AmBe_TH80_HV1550_D1md5cm_400ns_24hrs/UNFILTERED/Data_CH0.bin";
const string ATriggerFileName = "20211226AmBe_TH80_HV1550_D1md5cm_400ns_24hrs/UNFILTERED/Data_CH4.bin";
const string BTriggerFileName = "20211226AmBe_TH80_HV1550_D1md5cm_400ns_24hrs/UNFILTERED/Data_CH6.bin";

const int Recordlength = 992;
const int bits = 192;
const int numberOfBytesOfEvent = 312;       // bits devided by 8 (8bits/byte)

int64_t neutron_count =0 ;


//Variable

int64_t fATrigNumofRead;
int64_t fBTrigNumofRead;
int64_t main_detector_signal_counter;
int64_t A_detector_signal_counter;
int64_t B_detector_signal_counter;

//int waveCounter;
//const int numOfWavesToLookAt = 100;

//Number of threads
int world_rank = 0;


//ATrigger EJ-309
ifstream ATriginFile;
ifstream fATrigIn;
vector<double> fATrigcurrentWave;
int16_t fATrigboard;
int16_t fATrigchannel;
int64_t fATrigtimestamp;
int16_t fATrigenergy;
int16_t fATrigshortEnergy;
int32_t fATrigflags;
int16_t fATrignumOfWaveSamplesToBeRead;
int16_t Asample;
double  fATrigPSP;
bool Aneutron = false;
bool AAneutron = false;
//int AEn;

//BTrigger EJ-309
ifstream BTriginFile;
ifstream fBTrigIn;
vector<double> fBTrigcurrentWave;
int16_t fBTrigboard;
int16_t fBTrigchannel;
int64_t fBTrigtimestamp;
int16_t fBTrigenergy;
int16_t fBTrigshortEnergy;
int32_t fBTrigflags;
int16_t fBTrignumOfWaveSamplesToBeRead;
int16_t Bsample;
double fBTrigPSP;
bool Bneutron = false;
bool BBneutron = false;
//int BEn;

//Main EJ-290
ifstream ROIinFile;
ifstream fROIIn;
ifstream MatrixinFile;


int16_t fROIboard;
int16_t fROIchannel;
int64_t fROItimestamp;
int16_t fROIenergy;
int16_t fROIshortEnergy;
int32_t fROIflags;
int16_t fROInumOfWaveSamplesToBeRead; //record length in samples
int16_t sample;
double  fROIlistPsp;
vector<double> fROIcurrentWave; //vector of current wave
double fROIlongGate, fROIshortGate, fROIbaseline, fROIpsp;
//int64_t fROITOF;

//File output

ofstream TrigoutFile, Trigwave;
ofstream BTrigoutFile, BTrigwave;
ofstream ROIoutFile, ROIwave;
ofstream BROIoutFile;
ofstream MaitrxoutFile;

string OutTrigFileName, OutTrigWaveFileName;
string OutBTrigFileName, OutBTrigWaveFileName;
string OutROIFileName, OutROIFileWaveName, OutMatrixFileName;
string OutBROIFileName;

//boundary
vector<double> g1u;
vector<double> g1d;
vector<double> g2u;
vector<double> g2d;
vector<double> nu;
vector<double> nd;

//Function
void ROIopenBinaryFile(string inputFileName);
void ROIreadBinaryWave(int64_t start_pos, int64_t end_pos, int world_rank);
void ROIcloseFile();
void WriteTo309File();
void WriteTo290File();

void ATrigopenBinaryFile(string inputFileName);
void ATrigreadBinaryWave();
void ATrigcloseFile();

void BTrigopenBinaryFile(string inputFileName);
void BTrigreadBinaryWave();
void BTrigcloseFile();

//TH1F *TrigwaveHisto[numOfWavesToLookAt];
//TH1F *MainwaveHisto[numOfWavesToLookAt];

void ROOTcreateWaveHistos();
void ROOTwriteToFile();
int64_t numberOfEvent;

double fROIPSP;

int main(int argc, char** argv){
	
	//read PSD fiducial cut boundary

	// read the gamma detector 1 fiducial cut boundary
	/*
    std::ifstream ifs1("g1_PSDboundary_1std.txt", std::ios::in);
    if (!ifs1.is_open()) {
        cout << "Failed to open file.\n";
    } else {
        float g1channel;
        float g1up;
		float g1down;
        while (ifs1 >> g1channel >> g1up >> g1down) {
            //cout << "g1" << g1channel<<" "<<g1up << " " << g1down << "\n";
            g1u.push_back(g1up);
            g1d.push_back(g1down);
			
        }
		
        ifs1.close();
    }
    ifs1.close();


	// read the gamma detector 2 fiducial cut boundary

    std::ifstream ifs2("g2_PSDboundary_1std.txt", std::ios::in);
    if (!ifs2.is_open()) {
        cout << "Failed to open file.\n";
    } else {
        float g2channel;
        float g2up;
		float g2down;
        while (ifs2 >> g2channel >> g2up >> g2down) {
            //cout << "g2" << g2channel <<" "<<g2up << " " << g2down << "\n";
            g2u.push_back(g2up);
            g2d.push_back(g2down);
			
        }
		
        ifs2.close();
    }
    ifs2.close();

	// read the neutron detector fiducial cut boundary

    std::ifstream ifs3("n_PSDboundary_1std.txt", std::ios::in);
    if (!ifs3.is_open()) {
        cout << "Failed to open file.\n";
    } else {
        float nchannel;
        float nup;
		float ndown;
        while (ifs3 >> nchannel >> nup >> ndown) {
            //cout << "n" << nchannel <<" "<<nup << " " << ndown << "\n";
            nu.push_back(nup);
            nd.push_back(ndown);
			
        }
		
        ifs3.close();
    }
    ifs3.close();
    
	cout << "done boubdary" << "\n";
	*/

	// read the file
	time_t time_start, time_end;
	time_start = time(NULL);

	ROIopenBinaryFile(MainFileName);
	ATrigopenBinaryFile(ATriggerFileName);
	BTrigopenBinaryFile(BTriggerFileName);
	
/*
	
	ROIinFile.seekg(0, ios::end);

	//MPI_Init(NULL, NULL);
	int world_size = 1;
	//MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	//int world_rank;
	//MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	numberOfEvent = ROIinFile.tellg() / numberOfBytesOfEvent;
	//int64_t numberOfEvent = 21001656;
	int64_t numberOfEventForEachThread = numberOfEvent / world_size;

	int64_t start_pos = world_rank * numberOfEventForEachThread * numberOfBytesOfEvent;
	//int64_t start_pos = 27882695250;
	int64_t end_pos = (world_rank + 1) * numberOfEventForEachThread * numberOfBytesOfEvent;
	cout << end_pos/numberOfBytesOfEvent << endl;

	if (world_rank == 0){
		cout << "Data size: " << ROIinFile.tellg() << endl;
		cout << "Number of events: " << numberOfEvent << endl;
		cout << "Number of processes: " << world_size << endl;
		cout << "Number of events for each processes: " << numberOfEventForEachThread << endl;
		cout << "Press enter to start" << endl;
		cin.get();
		//MPI_Bcast(&TrigtotalEventCounter, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else{
		//MPI_Bcast(&TrigtotalEventCounter, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

	stringstream ss;
	ss << world_rank;
	string NumofFile = ss.str();

	OutROIFileName =  "EJ290@0-250_EN" + NumofFile + ".txt";
	OutBROIFileName = "EJ290@50-200_EN" + NumofFile + ".txt";
	//OutTrigFileName = "EJ309_-80ns" + NumofFile + ".txt";
	//OutBTrigFileName = "EJ309_+200ns" + NumofFile + ".txt";
	
	//OutROIFileWaveName = "aMainSignalWave" + NumofFile + ".txt";
	//OutTrigWaveFileName = "aTrigWave" + NumofFile + ".txt";
	//OutVeriWaveFileName = "aVeriWave" + NumofFile + ".txt";
	
	ROIoutFile.open(OutROIFileName, std::ios_base::app);
	BROIoutFile.open(OutBROIFileName,std::ios_base::app);
	//TrigoutFile.open(OutTrigFileName, std::ios_base::app);
	//BTrigoutFile.open(OutBTrigFileName, std::ios_base::app);
	
	//ROIwave.open(OutROIFileWaveName, std::ios_base::app);
	//Trigwave.open(OutTrigWaveFileName, std::ios_base::app);
	//Veriwave.open(OutVeriWaveFileName, std::ios_base::app);
	
	//waveCounter = 0 ;
*/

	// open output file for write //
	stringstream ss;
	ss << world_rank;
	string NumofFile = ss.str();

	OutROIFileName =  "20211227_TOF_D1m_d5cm_N_woPSD_g1g2_wPSD_24hrs.txt";  //need change
	/*
	Rule of output file name
	D 3" EJ309 to source
	d 2" EJ309 to source
	N: neutron ; g1,g2 gamma trigger detector
	wo,w : without, with
	*/
	OutMatrixFileName = "20211227_time_overlay_24hrs.txt";
	//OutROIFileName =  "EJ290@70-140_EN" + NumofFile + ".txt";
	//OutBROIFileName = "EJ290@50-200_EN" + NumofFile + ".txt";
	
	ROIoutFile.open(OutROIFileName, std::ios_base::app);
	MaitrxoutFile.open(OutMatrixFileName, std::ios_base::app);
	//BROIoutFile.open(OutBROIFileName,std::ios_base::app);
	
	
	//ROOTcreateWaveHistos();
	int64_t start_pos = 0;
	int64_t end_pos = 0;
	ROIreadBinaryWave(start_pos, end_pos, world_rank);


	//ROIinFile.seekg(0);

	ATrigcloseFile();
	BTrigcloseFile();
	ROIcloseFile();

	//cout<<"End of analysis"<<endl;


	MaitrxoutFile  <<  "test" << endl;

	time_end = time(NULL);
	cout << time_end - time_start << " second" << endl;

	cout << neutron_count << " neutrons" << endl;
	
	//ROOTwriteToFile();
	

	// MPI_Finalize();
	return 0;
}

/*
void WriteToFile(){

	if ((Aneutron == true) && (Bneutron == false)){

		neutron_count++;
		//cout << "ATRIG" << endl;
		
		int64_t ROITIME = fROItimestamp/1000/2;
		int64_t ATRIGTIME = fATrigtimestamp/1000/2;
		int TIME = 380+ATRIGTIME-ROITIME; 


		for(int i=TIME ; i < fATrigcurrentWave.size()+TIME ; i++){ 
			
    		TrigwaveHisto[waveCounter]->SetBinContent(i,fATrigcurrentWave.at(i-TIME));
    	}
		for(int i= 380 ; i < fROIcurrentWave.size()+380  ; i++){ 
    		MainwaveHisto[waveCounter]->SetBinContent(i,fROIcurrentWave.at(i-380));
    	}
    	waveCounter++;
 	}
	if ((Aneutron == false) && (Bneutron == true)){
		neutron_count++;
		//cout << "BTRIG" << endl;
		
		int64_t ROITIME = fROItimestamp/1000/2;
		int64_t BTRIGTIME = fBTrigtimestamp/1000/2;
		int TIME = 380+BTRIGTIME-ROITIME; 

		for(int i=TIME; i< fBTrigcurrentWave.size()+TIME; i++){ 

			//cout << "BTRIG11" << endl;
			//cout << fBTrigcurrentWave.size() << endl;
			//cout << waveCounter << endl;

    		TrigwaveHisto[waveCounter]->SetBinContent(i,fBTrigcurrentWave.at(i-TIME));

			//cout << i << endl;
    	}
		//cout << "ENDOFTRIG" << endl;
		
		for(int i=380; i< fROIcurrentWave.size()+380 ;i++){ 

			//cout << "BTRIG22" << endl;
    		MainwaveHisto[waveCounter]->SetBinContent(i,fROIcurrentWave.at(i-380));

    	}
		
    	waveCounter++;
 	}
	// cout << waveCounter<< endl;
	
}

*/
void ROOTwriteToFile(){
	/*
	TFile *run = new TFile("waves.root","recreate");
	run->SetCompressionSettings(3);

  	for(int i=0; i<numOfWavesToLookAt; i++)
  	{
    	//waveHisto[i]->SetLineColor(i);
    	TCanvas *c1 = new TCanvas("c1");

		MainwaveHisto[i]->SetLineColor(1);
		MainwaveHisto[i]->Draw();

		TrigwaveHisto[i]->SetLineColor(2);
		TrigwaveHisto[i]->Draw("same");


		
		c1->Write();
		delete c1;
  	}
	
  	run->Close();  
*/	
}

void ROOTcreateWaveHistos(){
	/*
	for(int i=0; i<numOfWavesToLookAt; i++){
    	stringstream SSname; //cout << i << endl;
    	TString histoName;// = SSname.str();
    	SSname << "wave_" << i;
    	histoName = SSname.str(); 
    	// 2 ns sampling is achieved by deselecting all other traces besides the 
    	// input in the waveform visualizatio gui
    	TrigwaveHisto[i] = new TH1F(histoName,";sample (2 ns);channel",800,-400,400);
		MainwaveHisto[i] = new TH1F(histoName,";sample (2 ns);channel",800,-400,400);
  }
  */
}

/////////////////EJ309-neutron//////////////////////////

void ROIopenBinaryFile(string inputFileName){
	//cout << inputFileName << endl;
	ROIinFile.open(inputFileName.c_str(), ios::binary | ios::in);
	if (!ROIinFile.good()) cout << "MainEJ290 file is not open" << '\n'; //check file is open

	//initize some things when you open the file
}



void ROIcloseFile(){
	ROIinFile.close();
}

void ROIreadBinaryWave(int64_t start_pos, int64_t end_pos, int nThread){
	

	
	fROIIn.open(MainFileName.c_str(), ios::binary | ios::in);
	fROIIn.seekg(start_pos, ios::beg);
	fATrigIn.open(ATriggerFileName.c_str(), ios::binary | ios::in);
	fBTrigIn.open(BTriggerFileName.c_str(), ios::binary | ios::in);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
	//fROIcurrentWave.clear();

	while (1){
		
		//Aneutron = false;
		
		/*
		if (waveCounter==100){
			break;
		}
		*/
		Aneutron = false;
		AAneutron = false;
		Bneutron = false;
		BBneutron = false;
		fROIcurrentWave.clear();
		fROIIn.read((char*)&fROIboard, sizeof(int16_t));
		fROIIn.read((char*)&fROIchannel, sizeof(int16_t));
		fROIIn.read((char*)&fROItimestamp, sizeof(int64_t)); //pico seconds
		fROIIn.read((char*)&fROIenergy, sizeof(int16_t));
		fROIIn.read((char*)&fROIshortEnergy, sizeof(int16_t));
		fROIIn.read((char*)&fROIflags, sizeof(int32_t));
		fROIIn.read((char*)&fROInumOfWaveSamplesToBeRead, sizeof(int16_t));
		fROIIn.read((char*)&sample, sizeof(int16_t));
		int16_t sample;
		for (int i = 0; i < fROInumOfWaveSamplesToBeRead; i++)
		{
			fROIIn.read((char*)&sample, sizeof(int16_t));
			fROIcurrentWave.push_back(sample);
		}
		
		
		
		//double fROIPSP;
		
		if (fROIenergy == 0) fROIPSP = 0;

		else { fROIPSP = (fROIenergy - fROIshortEnergy) / ((double)fROIenergy *1.0); }


		
		// display term 
		/*
		cout << "__________________mainEJ290______________" << endl;
		cout << "This is for the first event of the file" << endl;
		cout << "board: " << fROIboard << endl;
		cout << "ch: " << fROIchannel << endl;
		cout << "ts: " << fROItimestamp << endl;
		cout << "E: " << fROIenergy << endl;
		cout << "sE: " << fROIshortEnergy << endl;
		cout << "PSP:" << fROIPSP << endl;
		cout << "rec. length: " << fROInumOfWaveSamplesToBeRead << endl;
		cout << "_________________________________________" << endl << endl;
		fgetc(stdin);
		*/

		// For fiducial cut
		/*	
		if (( fROIPSP > (nd[fROIenergy-1]/1000) )&&( fROIPSP < (nu[fROIenergy-1]/1000))){
			neutron_count++;
		//if (( fROIPSP > 0.31 ) && ( fROIPSP < 0.44 )){
		//if ( fROIPSP > 0 ) {
			//cout << "okayyyyyyy"<< endl;
			ATrigreadBinaryWave();
			
		}
		*/
		
		if ( fROIPSP < 0 ){
			continue;
		}


		//if (( fROIPSP > 0.3 ) && ( fROIPSP < 0.5 )){
		neutron_count++;
		ATrigreadBinaryWave();
		
		//}


		/*
		bool checkFile;
		checkFile = fROIIn.eof(); //this will be true at the end of the file
		if (checkFile){
			cout << " B finished "<<endl;
			return;
		}
		*/
		/*
		ATrigreadBinaryWave();
		*/

		
		
		main_detector_signal_counter++;
		if (main_detector_signal_counter % 10000 == 0){
			cout << "Event: "<<main_detector_signal_counter << endl<<endl;;
		
		}
		/*
		if (neutron_count % 10000 == 0){
			cout << "Neutron counts: "<<neutron_count << endl<<endl;;
		
		}
		*/
		/*
		if (main_detector_signal_counter == 10000000){
			cout<<"end of the run at "<< main_detector_signal_counter <<endl;
			break;
		}
		*/
		
		/*
		if (main_detector_signal_counter >= 6283000){ //
			cout << main_detector_signal_counter << endl;
			cout << "time "<< fROItimestamp << endl;
			fgetc(stdin);
			//return;
		}
		
		
		if (fROItimestamp >= 288000000000000000){ //6283068
			cout << main_detector_signal_counter << endl;
			cout << "time "<< fROItimestamp << endl;
			fgetc(stdin);
			//return;
		}
		*/


		if(!fROIIn){
			cout<<"break in main loop"<< endl;
			cout << "Neutron counts: "<<neutron_count << endl<<endl;;
			cout << "Event: "<<main_detector_signal_counter << endl<<endl;;
			fROIIn.close();
			fATrigIn.close();
			fBTrigIn.close();
			break;
		}		

		/*
		if (fROIIn.tellg() >= end_pos || fROIIn.eof()){
			fROIIn.close();
			fATrigIn.close();
			fBTrigIn.close();
			cout << "End of thread: " << nThread << endl;
			return;
		}
		*/
		
		}
		if(!fROIIn){
		cout<<"break in outside loop"<< endl;
		cout << "Neutron counts: "<<neutron_count << endl<<endl;;
		cout << "Event: "<<main_detector_signal_counter << endl<<endl;;
		fROIIn.close();
		fATrigIn.close();
		fBTrigIn.close();
	}		
	return;
}

/////////////////EJ309-1//////////////////////////

void ATrigopenBinaryFile(string inputFileName){
	//cout << inputFileName << endl;
	ATriginFile.open(inputFileName.c_str(), ios::binary | ios::in);
	if (!ATriginFile.good()) cout << "Trigger A file is not open" << '\n'; //check file is open

	//initize some things when you open the file
	//TrigtotalEventCounter = 0;
	//TrigwaveCounter = 0;

}


void ATrigcloseFile(){
	ATriginFile.close();
}


void ATrigreadBinaryWave(){

	//fATrigcurrentWave.clear();
	int64_t positionA = numberOfBytesOfEvent* fATrigNumofRead;
	fATrigIn.seekg(positionA, ios::beg);

	while (1){
		
		//Aneutron = false;
		//AEn = -1;
		Aneutron = false;
		AAneutron = false;
		fATrigcurrentWave.clear();
		fATrigIn.read((char*)&fATrigboard, sizeof(int16_t));
		fATrigIn.read((char*)&fATrigchannel, sizeof(int16_t));
		fATrigIn.read((char*)&fATrigtimestamp, sizeof(int64_t)); //pico seconds
		fATrigIn.read((char*)&fATrigenergy, sizeof(int16_t));
		fATrigIn.read((char*)&fATrigshortEnergy, sizeof(int16_t));
		fATrigIn.read((char*)&fATrigflags, sizeof(int32_t));
		fATrigIn.read((char*)&fATrignumOfWaveSamplesToBeRead, sizeof(int16_t));
		fATrigIn.read((char*)&Asample, sizeof(int16_t));

	
		
		// read waveform
		
		int16_t Asample;
		for (int i = 0; i < fATrignumOfWaveSamplesToBeRead; i++)
		{
			fATrigIn.read((char*)&Asample, sizeof(int16_t));
			fATrigcurrentWave.push_back(Asample);
		}

		//cout << "fine"<< endl;
		if (fATrigenergy == 0) fATrigPSP = 0;
		
		else { fATrigPSP = (fATrigenergy - fATrigshortEnergy) / ((double)fATrigenergy *1.0); }
		//else { fATrigPSP = 1-(fATrigshortEnergy/(double)fATrigenergy); }
		/*
		cout << "_________________TriggerA______________" << endl;
		cout << "This is for the first event of the file" << endl;
		cout << "board: " << fATrigboard << endl;
		cout << "ch: " << fATrigchannel << endl;
		cout << "ts: " << fATrigtimestamp << endl;
		cout << "E: " << fATrigenergy << endl;
		cout << "sE: " << fATrigshortEnergy << endl;
		cout << "Flag " << fATrigflags << endl;
		cout << "rec. length: " << fATrignumOfWaveSamplesToBeRead << endl;
		cout << "PSP: " << fATrigPSP << endl;
		cout << "DeltaT: " << (fROItimestamp - fATrigtimestamp)/1000 << endl;
		cout << "T290: " << fROItimestamp / 1000 << endl;
		cout << "fATrigNumofRead: " << fATrigNumofRead << endl;
		cout << "_________________________________________" << endl;

		fgetc(stdin); //pause sign
		*/

		double TIMEDIV = ((fROItimestamp*1.0) - (fATrigtimestamp*1.0) ) / (1000.0) ;
		
			//cout << "A "<<TIMEDIV << endl;
			//cout << "APSP "<<fATrigPSP<<endl;
			//fgetc(stdin);
			
		///// if ((fATrigPSP > (g1d[fATrigenergy-1]/1000)) && (fATrigPSP < (g1u[fATrigenergy-1]/1000)) )  {
		
		
		//if ((fATrigPSP > 0.1) && (fATrigPSP < 0.28))  {
		//if (fATrigPSP > 0)   {
			
			// determine the available
			/*
			if ((TIMEDIV<=350)&&(TIMEDIV>=-150)){
				AAneutron = true;
				BTrigreadBinaryWave();
				fATrigNumofRead++;
				break;
				
			}
		
			if (TIMEDIV<-400){
				BTrigreadBinaryWave();
				break;
			}
			if (TIMEDIV>400){
				//BTrigreadBinaryWave();
				fATrigNumofRead++;
				cout <<'time diviation larger than positive 400'<<endl;
			}
			else {
				//cout << "WRONG AGAIN"<< endl;
				BTrigreadBinaryWave();
				fATrigNumofRead++;
				//break;
			}
			*/
		// if ((fATrigPSP > (g1d[fATrigenergy-1]/1000)) && (fATrigPSP < (g1u[fATrigenergy-1]/1000)) )  {
		
		if ((fATrigPSP > 0.1) && (fATrigPSP < 0.3))  {

		if ((TIMEDIV<=399)&&(TIMEDIV>=-399)){
			AAneutron = true;
			BTrigreadBinaryWave();
			fATrigNumofRead++;
			break;
				
			}
			
		if (TIMEDIV<-400){
			BTrigreadBinaryWave();
			break;
			}
		if (TIMEDIV>400){
			//BTrigreadBinaryWave();
			fATrigNumofRead++;
			}
		else {
			//cout << "WRONG AGAIN"<< endl;
			BTrigreadBinaryWave();
			fATrigNumofRead++;
			//break;
			}
		}

		/*
		if (TIMEDIV<-400){
			BTrigreadBinaryWave();
			//cout <<"A  "<<TIMEDIV<< "   time diviation lower than positive -400"<< endl;
			break;
		}
		if (TIMEDIV>400){
			//BTrigreadBinaryWave();
			fATrigNumofRead++;
			//cout <<"A  "<<TIMEDIV<< "   time diviation larger than positive 400"<< endl;
		}
		else {
			//cout << "WRONG AGAIN"<< endl;
			AAneutron = true;
			BTrigreadBinaryWave();
			fATrigNumofRead++;
			break;
			//cout <<"A  "<<TIMEDIV<< "  triggered  "<< endl;
			//fgetc(stdin); //pause sign
			//break;
		}
		*/
		///// }
		
		///// else{
		///// 	//BTrigreadBinaryWave();
		///// 	fATrigNumofRead++;
		///// 	//break;
		///// }
		
		
		//fgetc(stdin);	
		/*
		if ((fATrigPSP <= 0.25) && (fATrigPSP >= 0)) {
			if ((TIMEDIV >= 50) && (TIMEDIV <= 400)){
				if ((TIMEDIV<=85)&&(TIMEDIV>=75)){
					AEn = 0;
				}
					Aneutron = true;
					fATrigNumofRead++;
					BTrigreadBinaryWave();
					return;
			}
			if (TIMEDIV < 50){
				if ((TIMEDIV<=-195)&&(TIMEDIV>=-205)){
					AEn = 1;
					//cout << "XXXXX"<<endl;
					//fgetc(stdin);
				}
				//cout << "22222" << endl;
				BTrigreadBinaryWave();
				return;
			}
			if (TIMEDIV > 400){
				fATrigNumofRead++;
				//cout << "33333" << endl;
				BTrigreadBinaryWave();
			}
			
			
		}
		else { // PSP>0.25 || PSP < 0 ,neutron || meaningless
			fATrigNumofRead ++;
			BTrigreadBinaryWave();
		}
		*/
		/*
		if(!fATrigIn){
			cout << " A finished inside"<<endl;
			break;
		}
		*/	


		/*
		bool checkFile;
		checkFile = fATrigIn.eof(); //this will be true at the end of the file
		if (checkFile){
			cout << " A finished "<<endl;
			return;
		}
		*/
		/*
		if(!fATrigIn){
		cout << " A finished at inside"<<endl;

		break;

		}
		*/
	}

	

	//cout << " out of reading A "<<endl;
	//fgetc(stdin);
	return;
}



/////////////////EJ309-2//////////////////////////

void BTrigopenBinaryFile(string inputFileName){
	//cout << inputFileName << endl;
	BTriginFile.open(inputFileName.c_str(), ios::binary | ios::in);
	if (!BTriginFile.good()) cout << "Trigger B file is not open" << '\n'; //check file is open

	//initize some things when you open the file
	//TrigtotalEventCounter = 0;
	//TrigwaveCounter = 0;

}

void BTrigcloseFile(){
	BTriginFile.close();
}


void BTrigreadBinaryWave(){

	fBTrigcurrentWave.clear();
	int64_t positionB = numberOfBytesOfEvent* fBTrigNumofRead;
	fBTrigIn.seekg(positionB, ios::beg);

	while (1){
		//Bneutron = false;
		//BEn = -1;
		Bneutron = false;
		BBneutron = false;
		fBTrigcurrentWave.clear();
		
		fBTrigIn.read((char*)&fBTrigboard, sizeof(int16_t));
		fBTrigIn.read((char*)&fBTrigchannel, sizeof(int16_t));
		fBTrigIn.read((char*)&fBTrigtimestamp, sizeof(int64_t)); //pico seconds
		fBTrigIn.read((char*)&fBTrigenergy, sizeof(int16_t));
		fBTrigIn.read((char*)&fBTrigshortEnergy, sizeof(int16_t));
		fBTrigIn.read((char*)&fBTrigflags, sizeof(int32_t));
		fBTrigIn.read((char*)&fBTrignumOfWaveSamplesToBeRead, sizeof(int16_t));
		fBTrigIn.read((char*)&Bsample, sizeof(int16_t));

	
		
		
	/*
		int16_t Bsample;
		
		for (int i = 0; i < fBTrignumOfWaveSamplesToBeRead; i++)
		{
			fBTrigIn.read((char*)&Bsample, sizeof(int16_t));
			fBTrigcurrentWave.push_back(Bsample);
		}
	*/
		//cout << "WTF"<< endl;
		if (fBTrigenergy == 0) fBTrigPSP = 0;
		
		else { fBTrigPSP = (fBTrigenergy - fBTrigshortEnergy) / ((double)fBTrigenergy *1.0); }

		/*
		cout << "_________________TriggerB______________" << endl;
		cout << "This is for the first event of the file" << endl;
		cout << "board: " << fBTrigboard << endl;
		cout << "ch: " << fBTrigchannel << endl;
		cout << "ts: " << fBTrigtimestamp << endl;
		cout << "E: " << fBTrigenergy << endl;
		cout << "sE: " << fBTrigshortEnergy << endl;
		cout << "rec. length: " << fBTrignumOfWaveSamplesToBeRead << endl;
		cout << "PSP: " << fBTrigPSP << endl;
		cout << "DeltaT: " << (fROItimestamp - fBTrigtimestamp) / 1000 << endl;
		cout << "T290: " << fROItimestamp / 1000 << endl;
		cout << "fBTrigNumofRead: " << fBTrigNumofRead << endl;
		cout << "_________________________________________" << endl;


		fgetc(stdin); //pause sign
		*/
		
		double TIMEDIV =((fROItimestamp*1.0) - (fBTrigtimestamp*1.0) ) / (1000.0);
		//cout << "BBBBB " << TIMEDIV << endl;

		// //	

		///// if ((fBTrigPSP > (g2d[fATrigenergy-1]/1000)) && (fBTrigPSP < (g2u[fATrigenergy-1]/1000)) )  {
	
		//if ((fBTrigPSP > 0.1) && (fBTrigPSP < 0.28))   {
		//if (fBTrigPSP > 0)    {
			
		//if ((fBTrigPSP > (g2d[fATrigenergy-1]/1000)) && (fBTrigPSP < (g2u[fATrigenergy-1]/1000)) )  {
		if ((fBTrigPSP > 0.1) && (fBTrigPSP < 0.3)) {
			if ((TIMEDIV<=399)&&(TIMEDIV>=-399)){
				BBneutron = true;
				WriteTo290File();
				fBTrigNumofRead++;

				break;				
			}
			if (TIMEDIV<-400){
				WriteTo290File();
				break;
			}
			if (TIMEDIV>400){
				fBTrigNumofRead++;
			}
			else {
				fBTrigNumofRead++;
				WriteTo290File();
				break;
			}
		}
		/*
		if (TIMEDIV<-400){
			WriteTo290File();
			break;
		}
		if (TIMEDIV>400){
			fBTrigNumofRead++;
			//cout <<"B  "<<TIMEDIV<< "   time diviation larger than positive 400 "<< endl;
		}
		else {
			BBneutron = true;
			WriteTo290File();
			fBTrigNumofRead++;
			//cout <<"B  "<<TIMEDIV<< "   triggered  "<< endl;
			//fgetc(stdin); //pause sign
			break;
		}
		*/


		///// }
		///// else{
		/////	fBTrigNumofRead++;
		/////	//break;
		///// }



		// //
		
			//cout << "B "<<TIMEDIV << endl;
			//cout << "BPSP "<<fBTrigPSP<<endl;
			//fgetc(stdin);

		/*
		if ((fBTrigPSP <= 0.25) && (fBTrigPSP >= 0)) {
			if ((TIMEDIV >= 50) && (TIMEDIV <= 400)){
				if ((TIMEDIV<=85)&&(TIMEDIV>=75)){
					BEn = 0;
				}
					Bneutron = true;
					fBTrigNumofRead++;
					if ((Aneutron == true) || (Bneutron == true)){
						WriteToFile();
						return;
					}
			}
			if (TIMEDIV < 50){
				if ((TIMEDIV<=-195)&&(TIMEDIV>=-205)){
					BEn = 1;
				}
				if ((Aneutron == false)&& ( Bneutron == false)){
					WriteToFile();
				}
				if ((Aneutron == true) || (Bneutron == true)){
					WriteToFile();
				}
				return;
			}
			if (TIMEDIV > 400){
				fBTrigNumofRead++;
			}
		}
		else{ // PSP>0.25 || PSP < 0 ,neutron || meaningless
			fBTrigNumofRead++;
		}
	
		*/
	/*
		if(!fBTrigIn){
			cout << " B finished inside"<<endl;
			break;
		}
	*/      
		/*
		bool checkFile;
		checkFile = fBTrigIn.eof(); //this will be true at the end of the file
		if (checkFile){
			cout << " B finished "<<endl;
			return;
		}
		*/
	/*
		if(!fBTrigIn){
		cout << " B finished at inside"<<endl;
		break;
		}
	*/
	}
	/*
	cout << " out of reading B "<<endl;
	fgetc(stdin);
	*/
	return;

}
void WriteTo290File(){
	/*
	if ((AAneutron == true)&&(BBneutron==false)){
		ROIoutFile  <<  ((fROItimestamp*1.0) - (fATrigtimestamp*1.0 )) / 1000.0 << endl;
	}
	
	if ((AAneutron == false)&&(BBneutron==true)){
		ROIoutFile  <<  ((fROItimestamp*1.0) - (fBTrigtimestamp*1.0 )) / 1000.0 << endl;
	}
	*/
	if ((AAneutron == true)&&(BBneutron==false)){
		ROIoutFile  <<  ((fROItimestamp*1.0) - (fATrigtimestamp*1.0 )) / 1000.0 <<"  "<< fROIPSP <<"  " << fROIenergy << endl;
	}
	
	if ((AAneutron == false)&&(BBneutron==true)){
		ROIoutFile  <<  ((fROItimestamp*1.0) - (fBTrigtimestamp*1.0 )) / 1000.0 <<"  " << fROIPSP  <<"  "<< fROIenergy <<endl;
	}

	if ((AAneutron == true)&&(BBneutron==true)){
		MaitrxoutFile  <<((fROItimestamp*1.0) - (fATrigtimestamp*1.0 )) / 1000.0 <<"  "<<  ((fROItimestamp*1.0) - (fBTrigtimestamp*1.0 )) / 1000.0 <<"  " << fROIPSP  <<"  "<< fROIenergy <<endl;
	}
}


/*
void WriteTo309File(){

	if ((AAneutron == true)&&(BBneutron==false)){
		BROIoutFile << fROIenergy << "," << ((fROItimestamp - fATrigtimestamp )) / 1000 << endl;

	}
	if ((AAneutron == false)&&(BBneutron==true)){
		BROIoutFile << fROIenergy << "," <<  ((fROItimestamp - fBTrigtimestamp )) / 1000 << endl;
	}
	

}
*/