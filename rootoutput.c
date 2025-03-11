Part1
#ifndef JPSIMDATATYPE_H
#define JPSIMDATATYPE_H

#include "JPSimStruct.hh"
#include <map>
#include <vector>
#include "TROOT.h"
#include "TimeStamp.hh"

typedef std::vector<JPSimHit> JPSimHitList;
typedef std::map<G4int /*PmtId*/, JPSimHitList> JPSimHitListMap;
typedef std::vector<JPSimHitListMap> JPSimHitListMapWithIntvl;
//typedef std::map<JPSimInterval, JPSimHitListMap> JPSimHitListMapWithIntvl;

//typedef std::map<G4int /*PmtId*/, std::map<G4int /* Segment Id */, std::vector<JPSimHit> > > JPSimTempHitList;
typedef std::map<G4int /*Segment Id*/, JPSimHitListMap> JPSimTempHitBuffer;
typedef std::map<G4int /*Segment Id*/, std::vector<JPSimTruthTree_t> > JPSimTempTruthBuffer;
typedef std::map<G4int /*Segment Id*/, std::vector<JPSimTrack_t> > JPSimTempTrackBuffer;

typedef std::pair<TimeStamp, TimeStamp> JPSimInterval;
typedef std::vector<std::pair<TimeStamp, TimeStamp> > JPSimIntervalList;

typedef std::vector<UInt_t> JPSimWaveform;
typedef std::map<G4int /*PmtId*/, JPSimWaveform> JPSimWaveformList;





#endif //JPSIMDATATYPE_H


Part2
// This file contains the Monte-Carlo tree structure of JPSim output file.

#ifndef JPSimOutput_H
#define JPSimOutput_H 1

#include "TROOT.h"
#include "ReadoutType.h"
#include <map>
#include <vector>
#include <set>
#include <iostream>


// ********** This is the RunHeader tree **********

struct JPSimRunHeader_t {
	Int_t RunId;
	Int_t WindowSize;		// In ns
	Int_t TriggerPosition;
	Double_t PhotonFactor;
	Double_t DynamicRange;	// -DynamicRange to 0.  In V
	Int_t Bit;
	Double_t DCOffset;
};


// ********** This is the SimTruthTree **********
//
struct JPSimPrimaryParticle_t {
	Int_t TrackId;
	Int_t PdgId;
	Double_t px;
	Double_t py;
	Double_t pz;
	Double_t Ek;
    JPSimPrimaryParticle_t()
    {
        TrackId = 0;
        PdgId = 0;
        px = 0;
        py = 0;
        pz = 0;
        Ek = 0;
    }
};

struct JPSimPE_t {
	Int_t PMTId;
	Int_t segmentId;
	Int_t primaryParticleId;
	Double_t photonX;
	Double_t photonY;
	Double_t photonZ;
	Double_t dETime;	 // Zero point is vertex time (0). In ns.
	Double_t photonTime; // Zero point is vertex time (0). In ns.
	Double_t HitTime;	 // Zero point is vertex time (0). In ns. No TT and TTS.
	Double_t PulseTime;  // Zero point is vertex time (0). In ns.
	Int_t PESec;
	Int_t PENanoSec;
    double PESubNanoSec;
	double HitPosInWindow;
	Double_t Wavelength;
	Int_t PEType;
    double Charge;
	JPSimPE_t() {
		PMTId = 0;
		primaryParticleId = 0;
		photonX = 0;
		photonY = 0;
		photonZ = 0;
		dETime = 0;
		photonTime = 0;
		HitTime = 0;
		PulseTime = 0;
		PESec = 0;
		PENanoSec = 0;
		PESubNanoSec = 0;
		HitPosInWindow = 0;
		Wavelength = 0;
		PEType = 0;
        Charge = 0;
	}
};

struct JPSimStepPoint_t
{
	Int_t nProcessType;
	Int_t nProcessSubType;
	Double_t fX;
	Double_t fY;
	Double_t fZ;
	Double_t fPx;
	Double_t fPy;
	Double_t fPz;
	Double_t fEk;
	Double_t fdE;
	Double_t fTime;
	Int_t nTargetZ;
	Int_t nTargetA;
	std::vector<Int_t> nSecondaryPdgId;
	JPSimStepPoint_t() {nTargetZ=-1; nTargetA=-1;}
};

struct JPSimTrack_t
{
	Int_t nSegmentId;
	Int_t nParentTrackId;
	Int_t nTrackId;
	Int_t nPrimaryId;
	Int_t nPdgId;
    Bool_t bDetectedPhoton;
	std::vector<JPSimStepPoint_t> StepPoints;
	JPSimTrack_t() {}
};

struct JPSimTruthTree_t {
	Int_t RunId;		// The ID of Run
	Int_t SegmentId;	// The ID of Segment
	Int_t VertexId;
	Int_t VertexRadZ;
	Int_t VertexRadA;
	Int_t nParticle;
	
	Double_t x;
	Double_t y;
	Double_t z;
	Int_t Sec;
	Int_t NanoSec;
	std::vector<JPSimPrimaryParticle_t> PrimaryParticleList;
	Double_t EkMerged;

	std::vector<Double_t> dEList;	// Energy deposit in tagged volumes

	std::vector<Double_t> userdefinedA;
	std::vector<Double_t> userdefinedB;

	Int_t nFiredPMT;

	Int_t CPh;			// The number of C photons
	Int_t SPh;			// The number of S photons
	Int_t APh;			// The number of all photons
	Int_t CPE;			// The number of C p.e.
	Int_t SPE;
	Int_t APE;

	std::vector<JPSimTrack_t> trackList;

	JPSimTruthTree_t()
	{
		nFiredPMT = 0;
		EkMerged = 0;
		CPE = 0;
		SPE = 0;
		APE = 0;
		APh = 0;
		CPh = 0;
		SPh = 0;
	}

};

constexpr int MAXPAR = 100;

// ********** This is the SimTriggerTree **********

struct JPSimTriggerInfoTree_t {
	std::vector<JPSimPE_t> PEList;
	std::vector<JPSimTruthTree_t> truthList;
	Readout_t  readout;
	JPSimTriggerInfoTree_t() {}
};

#endif


Part3
#ifndef JPSIMSTRUCT_H
#define JPSIMSTRUCT_H

#include "G4RotationMatrix.hh"
#include "TimeStamp.hh"
#include <vector>
#include "G4SystemOfUnits.hh"
#include "TROOT.h"

struct TimeStampDouble
{
    TimeStamp t0;
    double t1;
};

struct JPSimHit
{
	int nSegmentId;	         // Simulate segment. not electronics segment
	int nParticleId;      
	//TimeStamp fHitTime;      // With TT and TTS. Global timestamp.
	TimeStampDouble PETime;      // With TT and TTS. Global timestamp.
	int nPmtId;
	int PEType;              // -1: Dark noise, 0: Cherekov, 1: Scin
	G4ThreeVector photonPos; // The photon emission position.
	double dETime;	         // Energy deposit time. Zero point is vertex time. In ns
	double photonTime;	     // Photon emission time.
	double HitTime;          // Without TT and TTS. Zero point is vertex time. In ns
	double PulseTime;        // With TT and TTS. Zero point is vertex time. In ns
	double Wavelength;       // In nm.
    double Charge;           // In ADC*ns.

	JPSimHit()
	{
		nSegmentId = -1;
		nParticleId = -1;
		//fHitTime = TimeStamp::GetBOT();
		PETime.t0 = TimeStamp::GetBOT();
        PETime.t1 = 0;
		nPmtId = -1;
		PEType = -1;
		photonPos = G4ThreeVector(0,0,0);
		dETime = -1;
		photonTime = -1;
		HitTime = -1;
		PulseTime = -1;
		Wavelength = -1;
        Charge = -1;
	}

	bool operator < (const JPSimHit& rhs) const { 
        double delta0 = this->PETime.t0-rhs.PETime.t0;
        double delta1 = this->PETime.t1-rhs.PETime.t1;
        return delta0+delta1<0;
        //return this->fHitTime < rhs.fHitTime; 
    } 

};

#endif //JPSIMSTRUCT_H


Part4
/*
  Output Tree structure
  Extracted from JPReadout.hh
  
  Zhe Wang, 2016, 10
*/

#ifndef _READOUT_TYPE_H_
#define _READOUT_TYPE_H_

struct Readout_t {
  Int_t RunNo;
  Int_t TriggerNo;
  Int_t TriggerType;
  Int_t DetectorID;
  Int_t Sec;
  Int_t NanoSec;
  std::vector<UInt_t> ChannelId;
  std::vector<UInt_t> Waveform;
};

#endif /* _READOUT_TYPE_H_ */






