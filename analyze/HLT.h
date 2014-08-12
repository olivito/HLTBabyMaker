// -*- C++ -*-
#ifndef HLT_H
#define HLT_H
#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"
#include <vector> 
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

using namespace std; 
class HLT {
private: 
protected: 
	unsigned int index;
	float calo_ht_;
	TBranch *calo_ht_branch;
	bool calo_ht_isLoaded;
	float calo_mht_eta_;
	TBranch *calo_mht_eta_branch;
	bool calo_mht_eta_isLoaded;
	float calo_mht_phi_;
	TBranch *calo_mht_phi_branch;
	bool calo_mht_phi_isLoaded;
	float calo_mht_pt_;
	TBranch *calo_mht_pt_branch;
	bool calo_mht_pt_isLoaded;
	float met_eta_;
	TBranch *met_eta_branch;
	bool met_eta_isLoaded;
	float met_phi_;
	TBranch *met_phi_branch;
	bool met_phi_isLoaded;
	float met_pt_;
	TBranch *met_pt_branch;
	bool met_pt_isLoaded;
	float pf_ht_;
	TBranch *pf_ht_branch;
	bool pf_ht_isLoaded;
	float pf_mht_eta_;
	TBranch *pf_mht_eta_branch;
	bool pf_mht_eta_isLoaded;
	float pf_mht_phi_;
	TBranch *pf_mht_phi_branch;
	bool pf_mht_phi_isLoaded;
	float pf_mht_pt_;
	TBranch *pf_mht_pt_branch;
	bool pf_mht_pt_isLoaded;
	float	scale1fb_;
	TBranch *scale1fb_branch;
	bool scale1fb_isLoaded;
	vector<float> calojets_eta_;
	TBranch *calojets_eta_branch;
	bool calojets_eta_isLoaded;
	vector<float> calojets_phi_;
	TBranch *calojets_phi_branch;
	bool calojets_phi_isLoaded;
	vector<float> calojets_pt_;
	TBranch *calojets_pt_branch;
	bool calojets_pt_isLoaded;
	vector<float> genjets_eta_;
	TBranch *genjets_eta_branch;
	bool genjets_eta_isLoaded;
	vector<float> genjets_phi_;
	TBranch *genjets_phi_branch;
	bool genjets_phi_isLoaded;
	vector<float> genjets_pt_;
	TBranch *genjets_pt_branch;
	bool genjets_pt_isLoaded;
	vector<float> pfjets_eta_;
	TBranch *pfjets_eta_branch;
	bool pfjets_eta_isLoaded;
	vector<float> pfjets_phi_;
	TBranch *pfjets_phi_branch;
	bool pfjets_phi_isLoaded;
	vector<float> pfjets_pt_;
	TBranch *pfjets_pt_branch;
	bool pfjets_pt_isLoaded;
public: 
void Init(TTree *tree) {
  tree->SetMakeClass(1);
	calo_ht_branch = 0;
	if (tree->GetAlias("calo_ht") != 0) {
		calo_ht_branch = tree->GetBranch(tree->GetAlias("calo_ht"));
		if (calo_ht_branch) {calo_ht_branch->SetAddress(&calo_ht_);}
	}
	calo_mht_eta_branch = 0;
	if (tree->GetAlias("calo_mht_eta") != 0) {
		calo_mht_eta_branch = tree->GetBranch(tree->GetAlias("calo_mht_eta"));
		if (calo_mht_eta_branch) {calo_mht_eta_branch->SetAddress(&calo_mht_eta_);}
	}
	calo_mht_phi_branch = 0;
	if (tree->GetAlias("calo_mht_phi") != 0) {
		calo_mht_phi_branch = tree->GetBranch(tree->GetAlias("calo_mht_phi"));
		if (calo_mht_phi_branch) {calo_mht_phi_branch->SetAddress(&calo_mht_phi_);}
	}
	calo_mht_pt_branch = 0;
	if (tree->GetAlias("calo_mht_pt") != 0) {
		calo_mht_pt_branch = tree->GetBranch(tree->GetAlias("calo_mht_pt"));
		if (calo_mht_pt_branch) {calo_mht_pt_branch->SetAddress(&calo_mht_pt_);}
	}
	met_eta_branch = 0;
	if (tree->GetAlias("met_eta") != 0) {
		met_eta_branch = tree->GetBranch(tree->GetAlias("met_eta"));
		if (met_eta_branch) {met_eta_branch->SetAddress(&met_eta_);}
	}
	met_phi_branch = 0;
	if (tree->GetAlias("met_phi") != 0) {
		met_phi_branch = tree->GetBranch(tree->GetAlias("met_phi"));
		if (met_phi_branch) {met_phi_branch->SetAddress(&met_phi_);}
	}
	met_pt_branch = 0;
	if (tree->GetAlias("met_pt") != 0) {
		met_pt_branch = tree->GetBranch(tree->GetAlias("met_pt"));
		if (met_pt_branch) {met_pt_branch->SetAddress(&met_pt_);}
	}
	pf_ht_branch = 0;
	if (tree->GetAlias("pf_ht") != 0) {
		pf_ht_branch = tree->GetBranch(tree->GetAlias("pf_ht"));
		if (pf_ht_branch) {pf_ht_branch->SetAddress(&pf_ht_);}
	}
	pf_mht_eta_branch = 0;
	if (tree->GetAlias("pf_mht_eta") != 0) {
		pf_mht_eta_branch = tree->GetBranch(tree->GetAlias("pf_mht_eta"));
		if (pf_mht_eta_branch) {pf_mht_eta_branch->SetAddress(&pf_mht_eta_);}
	}
	pf_mht_phi_branch = 0;
	if (tree->GetAlias("pf_mht_phi") != 0) {
		pf_mht_phi_branch = tree->GetBranch(tree->GetAlias("pf_mht_phi"));
		if (pf_mht_phi_branch) {pf_mht_phi_branch->SetAddress(&pf_mht_phi_);}
	}
	pf_mht_pt_branch = 0;
	if (tree->GetAlias("pf_mht_pt") != 0) {
		pf_mht_pt_branch = tree->GetBranch(tree->GetAlias("pf_mht_pt"));
		if (pf_mht_pt_branch) {pf_mht_pt_branch->SetAddress(&pf_mht_pt_);}
	}
	scale1fb_branch = 0;
	if (tree->GetAlias("scale1fb") != 0) {
		scale1fb_branch = tree->GetBranch(tree->GetAlias("scale1fb"));
		if (scale1fb_branch) {scale1fb_branch->SetAddress(&scale1fb_);}
	}
	calojets_eta_branch = 0;
	if (tree->GetAlias("calojets_eta") != 0) {
		calojets_eta_branch = tree->GetBranch(tree->GetAlias("calojets_eta"));
		if (calojets_eta_branch) {calojets_eta_branch->SetAddress(&calojets_eta_);}
	}
	calojets_phi_branch = 0;
	if (tree->GetAlias("calojets_phi") != 0) {
		calojets_phi_branch = tree->GetBranch(tree->GetAlias("calojets_phi"));
		if (calojets_phi_branch) {calojets_phi_branch->SetAddress(&calojets_phi_);}
	}
	calojets_pt_branch = 0;
	if (tree->GetAlias("calojets_pt") != 0) {
		calojets_pt_branch = tree->GetBranch(tree->GetAlias("calojets_pt"));
		if (calojets_pt_branch) {calojets_pt_branch->SetAddress(&calojets_pt_);}
	}
	genjets_eta_branch = 0;
	if (tree->GetAlias("genjets_eta") != 0) {
		genjets_eta_branch = tree->GetBranch(tree->GetAlias("genjets_eta"));
		if (genjets_eta_branch) {genjets_eta_branch->SetAddress(&genjets_eta_);}
	}
	genjets_phi_branch = 0;
	if (tree->GetAlias("genjets_phi") != 0) {
		genjets_phi_branch = tree->GetBranch(tree->GetAlias("genjets_phi"));
		if (genjets_phi_branch) {genjets_phi_branch->SetAddress(&genjets_phi_);}
	}
	genjets_pt_branch = 0;
	if (tree->GetAlias("genjets_pt") != 0) {
		genjets_pt_branch = tree->GetBranch(tree->GetAlias("genjets_pt"));
		if (genjets_pt_branch) {genjets_pt_branch->SetAddress(&genjets_pt_);}
	}
	pfjets_eta_branch = 0;
	if (tree->GetAlias("pfjets_eta") != 0) {
		pfjets_eta_branch = tree->GetBranch(tree->GetAlias("pfjets_eta"));
		if (pfjets_eta_branch) {pfjets_eta_branch->SetAddress(&pfjets_eta_);}
	}
	pfjets_phi_branch = 0;
	if (tree->GetAlias("pfjets_phi") != 0) {
		pfjets_phi_branch = tree->GetBranch(tree->GetAlias("pfjets_phi"));
		if (pfjets_phi_branch) {pfjets_phi_branch->SetAddress(&pfjets_phi_);}
	}
	pfjets_pt_branch = 0;
	if (tree->GetAlias("pfjets_pt") != 0) {
		pfjets_pt_branch = tree->GetBranch(tree->GetAlias("pfjets_pt"));
		if (pfjets_pt_branch) {pfjets_pt_branch->SetAddress(&pfjets_pt_);}
	}
  tree->SetMakeClass(0);
}
void GetEntry(unsigned int idx) 
	// this only marks branches as not loaded, saving a lot of time
	{
		index = idx;
		calo_ht_isLoaded = false;
		calo_mht_eta_isLoaded = false;
		calo_mht_phi_isLoaded = false;
		calo_mht_pt_isLoaded = false;
		met_eta_isLoaded = false;
		met_phi_isLoaded = false;
		met_pt_isLoaded = false;
		pf_ht_isLoaded = false;
		pf_mht_eta_isLoaded = false;
		pf_mht_phi_isLoaded = false;
		pf_mht_pt_isLoaded = false;
		scale1fb_isLoaded = false;
		calojets_eta_isLoaded = false;
		calojets_phi_isLoaded = false;
		calojets_pt_isLoaded = false;
		genjets_eta_isLoaded = false;
		genjets_phi_isLoaded = false;
		genjets_pt_isLoaded = false;
		pfjets_eta_isLoaded = false;
		pfjets_phi_isLoaded = false;
		pfjets_pt_isLoaded = false;
	}

void LoadAllBranches() 
	// load all branches
{
	if (calo_ht_branch != 0) calo_ht();
	if (calo_mht_eta_branch != 0) calo_mht_eta();
	if (calo_mht_phi_branch != 0) calo_mht_phi();
	if (calo_mht_pt_branch != 0) calo_mht_pt();
	if (met_eta_branch != 0) met_eta();
	if (met_phi_branch != 0) met_phi();
	if (met_pt_branch != 0) met_pt();
	if (pf_ht_branch != 0) pf_ht();
	if (pf_mht_eta_branch != 0) pf_mht_eta();
	if (pf_mht_phi_branch != 0) pf_mht_phi();
	if (pf_mht_pt_branch != 0) pf_mht_pt();
	if (scale1fb_branch != 0) scale1fb();
	if (calojets_eta_branch != 0) calojets_eta();
	if (calojets_phi_branch != 0) calojets_phi();
	if (calojets_pt_branch != 0) calojets_pt();
	if (genjets_eta_branch != 0) genjets_eta();
	if (genjets_phi_branch != 0) genjets_phi();
	if (genjets_pt_branch != 0) genjets_pt();
	if (pfjets_eta_branch != 0) pfjets_eta();
	if (pfjets_phi_branch != 0) pfjets_phi();
	if (pfjets_pt_branch != 0) pfjets_pt();
}

	float &calo_ht()
	{
		if (not calo_ht_isLoaded) {
			if (calo_ht_branch != 0) {
				calo_ht_branch->GetEntry(index);
			} else { 
				printf("branch calo_ht_branch does not exist!\n");
				exit(1);
			}
			calo_ht_isLoaded = true;
		}
		return calo_ht_;
	}
	float &calo_mht_eta()
	{
		if (not calo_mht_eta_isLoaded) {
			if (calo_mht_eta_branch != 0) {
				calo_mht_eta_branch->GetEntry(index);
			} else { 
				printf("branch calo_mht_eta_branch does not exist!\n");
				exit(1);
			}
			calo_mht_eta_isLoaded = true;
		}
		return calo_mht_eta_;
	}
	float &calo_mht_phi()
	{
		if (not calo_mht_phi_isLoaded) {
			if (calo_mht_phi_branch != 0) {
				calo_mht_phi_branch->GetEntry(index);
			} else { 
				printf("branch calo_mht_phi_branch does not exist!\n");
				exit(1);
			}
			calo_mht_phi_isLoaded = true;
		}
		return calo_mht_phi_;
	}
	float &calo_mht_pt()
	{
		if (not calo_mht_pt_isLoaded) {
			if (calo_mht_pt_branch != 0) {
				calo_mht_pt_branch->GetEntry(index);
			} else { 
				printf("branch calo_mht_pt_branch does not exist!\n");
				exit(1);
			}
			calo_mht_pt_isLoaded = true;
		}
		return calo_mht_pt_;
	}
	float &met_eta()
	{
		if (not met_eta_isLoaded) {
			if (met_eta_branch != 0) {
				met_eta_branch->GetEntry(index);
			} else { 
				printf("branch met_eta_branch does not exist!\n");
				exit(1);
			}
			met_eta_isLoaded = true;
		}
		return met_eta_;
	}
	float &met_phi()
	{
		if (not met_phi_isLoaded) {
			if (met_phi_branch != 0) {
				met_phi_branch->GetEntry(index);
			} else { 
				printf("branch met_phi_branch does not exist!\n");
				exit(1);
			}
			met_phi_isLoaded = true;
		}
		return met_phi_;
	}
	float &met_pt()
	{
		if (not met_pt_isLoaded) {
			if (met_pt_branch != 0) {
				met_pt_branch->GetEntry(index);
			} else { 
				printf("branch met_pt_branch does not exist!\n");
				exit(1);
			}
			met_pt_isLoaded = true;
		}
		return met_pt_;
	}
	float &pf_ht()
	{
		if (not pf_ht_isLoaded) {
			if (pf_ht_branch != 0) {
				pf_ht_branch->GetEntry(index);
			} else { 
				printf("branch pf_ht_branch does not exist!\n");
				exit(1);
			}
			pf_ht_isLoaded = true;
		}
		return pf_ht_;
	}
	float &pf_mht_eta()
	{
		if (not pf_mht_eta_isLoaded) {
			if (pf_mht_eta_branch != 0) {
				pf_mht_eta_branch->GetEntry(index);
			} else { 
				printf("branch pf_mht_eta_branch does not exist!\n");
				exit(1);
			}
			pf_mht_eta_isLoaded = true;
		}
		return pf_mht_eta_;
	}
	float &pf_mht_phi()
	{
		if (not pf_mht_phi_isLoaded) {
			if (pf_mht_phi_branch != 0) {
				pf_mht_phi_branch->GetEntry(index);
			} else { 
				printf("branch pf_mht_phi_branch does not exist!\n");
				exit(1);
			}
			pf_mht_phi_isLoaded = true;
		}
		return pf_mht_phi_;
	}
	float &pf_mht_pt()
	{
		if (not pf_mht_pt_isLoaded) {
			if (pf_mht_pt_branch != 0) {
				pf_mht_pt_branch->GetEntry(index);
			} else { 
				printf("branch pf_mht_pt_branch does not exist!\n");
				exit(1);
			}
			pf_mht_pt_isLoaded = true;
		}
		return pf_mht_pt_;
	}
	float &scale1fb()
	{
		if (not scale1fb_isLoaded) {
			if (scale1fb_branch != 0) {
				scale1fb_branch->GetEntry(index);
			} else { 
				printf("branch scale1fb_branch does not exist!\n");
				exit(1);
			}
			scale1fb_isLoaded = true;
		}
		return scale1fb_;
	}
	const vector<float> &calojets_eta()
	{
		if (not calojets_eta_isLoaded) {
			if (calojets_eta_branch != 0) {
				calojets_eta_branch->GetEntry(index);
			} else { 
				printf("branch calojets_eta_branch does not exist!\n");
				exit(1);
			}
			calojets_eta_isLoaded = true;
		}
		return calojets_eta_;
	}
	const vector<float> &calojets_phi()
	{
		if (not calojets_phi_isLoaded) {
			if (calojets_phi_branch != 0) {
				calojets_phi_branch->GetEntry(index);
			} else { 
				printf("branch calojets_phi_branch does not exist!\n");
				exit(1);
			}
			calojets_phi_isLoaded = true;
		}
		return calojets_phi_;
	}
	const vector<float> &calojets_pt()
	{
		if (not calojets_pt_isLoaded) {
			if (calojets_pt_branch != 0) {
				calojets_pt_branch->GetEntry(index);
			} else { 
				printf("branch calojets_pt_branch does not exist!\n");
				exit(1);
			}
			calojets_pt_isLoaded = true;
		}
		return calojets_pt_;
	}
	const vector<float> &genjets_eta()
	{
		if (not genjets_eta_isLoaded) {
			if (genjets_eta_branch != 0) {
				genjets_eta_branch->GetEntry(index);
			} else { 
				printf("branch genjets_eta_branch does not exist!\n");
				exit(1);
			}
			genjets_eta_isLoaded = true;
		}
		return genjets_eta_;
	}
	const vector<float> &genjets_phi()
	{
		if (not genjets_phi_isLoaded) {
			if (genjets_phi_branch != 0) {
				genjets_phi_branch->GetEntry(index);
			} else { 
				printf("branch genjets_phi_branch does not exist!\n");
				exit(1);
			}
			genjets_phi_isLoaded = true;
		}
		return genjets_phi_;
	}
	const vector<float> &genjets_pt()
	{
		if (not genjets_pt_isLoaded) {
			if (genjets_pt_branch != 0) {
				genjets_pt_branch->GetEntry(index);
			} else { 
				printf("branch genjets_pt_branch does not exist!\n");
				exit(1);
			}
			genjets_pt_isLoaded = true;
		}
		return genjets_pt_;
	}
	const vector<float> &pfjets_eta()
	{
		if (not pfjets_eta_isLoaded) {
			if (pfjets_eta_branch != 0) {
				pfjets_eta_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_eta_branch does not exist!\n");
				exit(1);
			}
			pfjets_eta_isLoaded = true;
		}
		return pfjets_eta_;
	}
	const vector<float> &pfjets_phi()
	{
		if (not pfjets_phi_isLoaded) {
			if (pfjets_phi_branch != 0) {
				pfjets_phi_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_phi_branch does not exist!\n");
				exit(1);
			}
			pfjets_phi_isLoaded = true;
		}
		return pfjets_phi_;
	}
	const vector<float> &pfjets_pt()
	{
		if (not pfjets_pt_isLoaded) {
			if (pfjets_pt_branch != 0) {
				pfjets_pt_branch->GetEntry(index);
			} else { 
				printf("branch pfjets_pt_branch does not exist!\n");
				exit(1);
			}
			pfjets_pt_isLoaded = true;
		}
		return pfjets_pt_;
	}

  static void progress( int nEventsTotal, int nEventsChain ){
    int period = 1000;
    if(nEventsTotal%1000 == 0) {
      // xterm magic from L. Vacavant and A. Cerri
      if (isatty(1)) {
        if( ( nEventsChain - nEventsTotal ) > period ){
          float frac = (float)nEventsTotal/(nEventsChain*0.01);
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
               "\033[0m\033[32m <---\033[0m\015", frac);
          fflush(stdout);
        }
        else {
          printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                 "\033[0m\033[32m <---\033[0m\015", 100.);
          cout << endl;
        }
      }
    }
  }
  
};

#ifndef __CINT__
extern HLT hlt;
#endif

namespace tas {
	const float &calo_ht();
	const float &calo_mht_eta();
	const float &calo_mht_phi();
	const float &calo_mht_pt();
	const float &met_eta();
	const float &met_phi();
	const float &met_pt();
	const float &pf_ht();
	const float &pf_mht_eta();
	const float &pf_mht_phi();
	const float &pf_mht_pt();
	const float &scale1fb();
	const vector<float> &calojets_eta();
	const vector<float> &calojets_phi();
	const vector<float> &calojets_pt();
	const vector<float> &genjets_eta();
	const vector<float> &genjets_phi();
	const vector<float> &genjets_pt();
	const vector<float> &pfjets_eta();
	const vector<float> &pfjets_phi();
	const vector<float> &pfjets_pt();
}
#endif
