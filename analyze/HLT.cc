#include "HLT.h"
HLT hlt;
namespace tas {
	const float &calo_ht() { return hlt.calo_ht(); }
	const float &calo_mht_eta() { return hlt.calo_mht_eta(); }
	const float &calo_mht_phi() { return hlt.calo_mht_phi(); }
	const float &calo_mht_pt() { return hlt.calo_mht_pt(); }
	const float &met_eta() { return hlt.met_eta(); }
	const float &met_phi() { return hlt.met_phi(); }
	const float &met_pt() { return hlt.met_pt(); }
	const float &pf_ht() { return hlt.pf_ht(); }
	const float &pf_mht_eta() { return hlt.pf_mht_eta(); }
	const float &pf_mht_phi() { return hlt.pf_mht_phi(); }
	const float &pf_mht_pt() { return hlt.pf_mht_pt(); }
	const float &scale1fb() { return hlt.scale1fb(); }
	const vector<float> &calojets_eta() { return hlt.calojets_eta(); }
	const vector<float> &calojets_phi() { return hlt.calojets_phi(); }
	const vector<float> &calojets_pt() { return hlt.calojets_pt(); }
	const vector<float> &genjets_eta() { return hlt.genjets_eta(); }
	const vector<float> &genjets_phi() { return hlt.genjets_phi(); }
	const vector<float> &genjets_pt() { return hlt.genjets_pt(); }
	const vector<float> &pfjets_eta() { return hlt.pfjets_eta(); }
	const vector<float> &pfjets_phi() { return hlt.pfjets_phi(); }
	const vector<float> &pfjets_pt() { return hlt.pfjets_pt(); }
}
