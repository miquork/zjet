// Same inversion, different tag policy and selected-probe cohort.
// Run only on productions containing ResponseAudit. This is a conditional
// four-flavor estimate, NOT a measurement of all data tagging transitions.
#include "analyzeFlavorMatrix.C"
#include "ZJetResponseAudit.h"

void compareTaggerResponse(const char *dataName="rootfiles/zjet_DATA.root",
                           const char *mcName="rootfiles/zjet_MC.root",
                           const char *outName="output/taggerRu/crossed_response.tsv",
                           bool mcClosure=false) {
  using namespace FlavorMatrixAnalysis;
  TFile data(dataName),mc(mcName);
  if(data.IsZombie()||mc.IsZombie()||!mc.Get("ResponseAudit/definition")||
     !data.Get("ResponseAudit/definition"))
    throw std::runtime_error("New ResponseAudit production required for crossed comparisons");
  std::ofstream out(outName);
  if(!out) throw std::runtime_error("Cannot create crossed-response table");
  out<<std::setprecision(12);
  out<<"cohort\tpolicy\tlow\thigh\tRu\tpurity_model\tflavor\tscale\terror_conditional\trank\tcondition\tunknown_mc_fraction\n";
  const double edges[]={30,50,70,100,150,210,300,400,600};
  const int ids[]={1,4,5,6}; const char *names[]={"uds","c","b","g"};
  const std::vector<std::vector<int>> truths={{1,3},{4},{5},{6}};
  for(const std::string cohort:{"new","legacy","common_new","common_legacy"})
    for(const auto &policy:ZJetResponseAudit::categories({})) {
      const std::string dir="ResponseAudit/"+cohort+"/"+policy.first+"/";
      auto *dc=dynamic_cast<TH3D*>(data.Get((dir+"counts").c_str()));
      auto *mcnt=dynamic_cast<TH3D*>(mc.Get((dir+"counts").c_str()));
      if(!dc||!mcnt) throw std::runtime_error("Missing crossed policy "+dir);
      for(int p=0;p<8;++p) {
        TMatrixD joint(4,4);TVectorD d(4),prior(4);joint.Zero();prior.Zero();
        double dt=0.,mt=0.,unknown=0.;
        for(int r=0;r<4;++r) {
          d[r]=mcClosure ? integratedCount(dc,ids[r],{1,3,4,5,6},edges[p],edges[p+1])
            : integratedDataCount(dc,ids[r],edges[p],edges[p+1]);dt+=d[r];
          unknown+=integratedCount(mcnt,ids[r],{0},edges[p],edges[p+1]);
          for(int f=0;f<4;++f) {
            const double n=integratedCount(mcnt,ids[r],truths[f],edges[p],edges[p+1]);
            if(n<0.) throw std::runtime_error("Negative joint probability; signed fit needed");
            joint(r,f)=n;prior[f]+=n;mt+=n;
          }
        }
        if(dt<=0.||mt<=0.) continue;
        joint*=1./mt;prior*=1./mt;d*=1./dt;
        for(const std::string model:{"mc","ipf"}) {
          TMatrixD inferred(joint);
          bool converged=model=="mc";
          if(model=="ipf") {
            for(int r=0;r<4;++r)for(int f=0;f<4;++f)inferred(r,f)+=1.e-12;
            for(int it=0;it<10000;++it) {
              for(int r=0;r<4;++r) {
                double total=0.;for(int f=0;f<4;++f)total+=inferred(r,f);
                if(total>0.)for(int f=0;f<4;++f)inferred(r,f)*=d[r]/total;
              }
              for(int f=0;f<4;++f) {
                double total=0.;for(int r=0;r<4;++r)total+=inferred(r,f);
                if(total>0.)for(int r=0;r<4;++r)inferred(r,f)*=prior[f]/total;
              }
              double error=0.;for(int r=0;r<4;++r) {
                double total=0.;for(int f=0;f<4;++f)total+=inferred(r,f);
                error=std::max(error,std::abs(total-d[r]));
              }
              if(error<1.e-10){converged=true;break;}
            }
          }
          if(!converged) continue;
          for(double ru:{.5,.92,1.1}) {
            TMatrixD a(4,4);TVectorD observed(4),errors(4);bool usable=true;
            for(int r=0;r<4;++r) {
              auto hdm=[&](TFile &file,const std::vector<int> &flavor) {
                ProfileSummary v[3];int i=0;
                for(const char *component:{"m0","mn","mu"})
                  v[i++]=summarizeProfile(dynamic_cast<TProfile3D*>(file.Get((dir+component).c_str())),
                                         ids[r],flavor,edges[p],edges[p+1]);
                ProfileSummary result;
                if(!v[0].valid||!v[1].valid||!v[2].valid)return result;
                const double den=1.-v[1].mean-v[2].mean/ru,num=v[0].mean-v[1].mean-v[2].mean;
                if(std::abs(den)<1.e-8)return result;
                result.mean=num/den;
                result.error=std::sqrt(std::pow(v[0].error/den,2)+
                  std::pow((num-den)*v[1].error/(den*den),2)+
                  std::pow((num/ru-den)*v[2].error/(den*den),2));
                result.valid=std::isfinite(result.mean);return result;
              };
              const auto obs=hdm(data,mcClosure ? std::vector<int>{1,3,4,5,6} : std::vector<int>{0});
              if(!obs.valid||obs.error<=0.){usable=false;break;}
              observed[r]=obs.mean;errors[r]=obs.error;
              // HDM(<components>) is nonlinear. Mix numerator and denominator,
              // not per-cell HDM with ordinary event purities. This design is
              // exact for flavor scales multiplying the cell numerator while
              // retaining MC cell mn,mu: an explicit, testable model assumption.
              double row=0.,den=0.;for(int f=0;f<4;++f)row+=inferred(r,f);
              if(row<=0.){usable=false;break;}
              for(int f=0;f<4;++f) {
                double values[3];int i=0;
                for(const char *component:{"m0","mn","mu"}) {
                  auto v=summarizeProfile(dynamic_cast<TProfile3D*>(mc.Get((dir+component).c_str())),
                                          ids[r],truths[f],edges[p],edges[p+1]);
                  values[i++]=v.valid?v.mean:NAN;
                }
                const double pur=inferred(r,f)/row;
                if(pur<1.e-10){a(r,f)=0.;continue;}
                if(!std::isfinite(values[0])||!std::isfinite(values[1])||!std::isfinite(values[2])){usable=false;break;}
                a(r,f)=pur*(values[0]-values[1]-values[2]);
                den+=pur*(1.-values[1]-values[2]/ru);
              }
              if(std::abs(den)<1.e-8){usable=false;break;}
              for(int f=0;f<4;++f)a(r,f)/=den;
            }
            if(!usable)continue;
            auto fit=fitResponse(a,observed,errors,0.); // No unity prior hiding instability.
            for(int f=0;f<4;++f)
              out<<cohort<<'\t'<<policy.first<<'\t'<<edges[p]<<'\t'<<edges[p+1]<<'\t'<<ru<<'\t'
                 <<model<<'\t'<<names[f]<<'\t'<<(fit.rank==4?fit.residual[f]:NAN)<<'\t'
                 <<(fit.rank==4?std::sqrt(std::max(0.,fit.covariance(f,f))):NAN)<<'\t'
                 <<fit.rank<<'\t'<<fit.fullCondition<<'\t'<<unknown/(unknown+mt)<<'\n';
          }
        }
      }
    }
  std::cout<<"Wrote conditional crossed fits to "<<outName<<std::endl;
}
