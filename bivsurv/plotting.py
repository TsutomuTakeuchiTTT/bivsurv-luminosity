"""Data-driven figures. Each file is one independent matplotlib figure.

No KDE, Gaussian smoothing, data-adaptive grid, or post-hoc CDF repair.
Maps display bin mass / dex^2. Grey/masked regions are NOT zero estimates.
Hatches show cells without two measured luminosities at the deepest limit.
"""
from __future__ import annotations
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Rectangle
from .model import GridModel
from .benchmark import bright_surface,save_json,write_csv

LABELS={'D_completion':'D-completion EM','no_D_ablation':'No-D ablation','A_only_histogram':'A-only histogram'}


def save_figure(fig,out,stem):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/f'{stem}.pdf',bbox_inches='tight',metadata={'Creator':'bivsurv-luminosity 0.1.0rc1','CreationDate':None,'ModDate':None})
    fig.savefig(out/f'{stem}.png',bbox_inches='tight',dpi=190)
    plt.close(fig)


def qmass(model,mass):
    x=mass.copy();x[model.invisible]=0
    return x/x.sum()


def make_all(output):
    root=Path(output);config=json.loads((root/'run_config.json').read_text());out=root/'figures'
    rows=list(csv.DictReader((root/'all_metrics.csv').open()))
    atlas=[];family_metadata={};representative_metrics=[]
    for rho in config['rho']:
        n=config['representative']['parent_size'];rep=config['representative']['replicate']
        case=f'rho{rho:.2f}_N{n}_r{rep:03d}'
        if not (root/'fits/deep'/f'{case}.npz').exists():continue
        for design in config['limits_standardized_X']:
            model=GridModel.from_dict(json.loads((root/'models'/f'{design}.json').read_text()))
            with np.load(root/'models'/f'{design}_arrays.npz') as v:
                x1,x2,we=[v[k].copy() for k in ('x_edges1','x_edges2','w_edges')]
            with np.load(root/'fits'/design/f'{case}.npz') as v:
                fit,noD,naive,truth=[v[k].copy() for k in ('mass','noD','naive','truth')]
            qfit,qtrue=qmass(model,fit),qmass(model,truth)
            cutoff=-model.limits_z[-1];prefix=f'rho{rho:.2f}_{design}'
            finite=(slice(1,-1),slice(1,-1))
            areas=np.diff(x1)[:,None]*np.diff(x2)[None,:]
            arrays={'truth_Q':qtrue,'naive_A':naive,'reconstructed_Q':qfit,'noD_Q':qmass(model,noD)}
            titles={'truth_Q':'True observable LF shape','naive_A':'Naive: both detections only','reconstructed_Q':'D-completion reconstruction','noD_Q':'Ablation: no truncation correction'}
            for tag,m in arrays.items():
                d=m.reshape(model.shape)[::-1,::-1][finite]/areas
                invisible=model.invisible.reshape(model.shape)[::-1,::-1][finite]
                resolved=model.common_detected.reshape(model.shape)[::-1,::-1][finite]
                image=np.ma.masked_where(~resolved | (d<=0),d)
                fig,ax=plt.subplots(figsize=(6.1,5.25))
                im=ax.pcolormesh(x1,x2,image.T,norm=LogNorm(vmin=1e-3,vmax=1.5),shading='flat',rasterized=True)
                cb=fig.colorbar(im,ax=ax,pad=.025);cb.set_label(r'Bin-averaged density [dex$^{-2}$]')
                ax.axvline(cutoff[0],linestyle=':',linewidth=1.1);ax.axhline(cutoff[1],linestyle=':',linewidth=1.1)
                # Unresolved strips are shown, not concealed; their lift is not point-identified.
                ax.add_patch(Rectangle((x1[0],x2[0]),cutoff[0]-x1[0],x2[-1]-x2[0],fill=False,hatch='/',linewidth=.25,alpha=.18))
                ax.add_patch(Rectangle((cutoff[0],x2[0]),x1[-1]-cutoff[0],cutoff[1]-x2[0],fill=False,hatch='/',linewidth=.25,alpha=.18))
                ax.text(x1[0]+.12, x2[0]+.22,'Invisible\nnot inferred',fontsize=8)
                ax.text(x1[0]+.12, cutoff[1]+.3,'Band 1 unresolved',fontsize=8,rotation=90)
                ax.text(cutoff[0]+.22, x2[0]+.18,'Band 2 unresolved',fontsize=8)
                ax.set(xlim=(x1[0],x1[-1]),ylim=(x2[0],x2[-1]),xlabel=r'$\log_{10}(L_1/L_{1,0})$',ylabel=r'$\log_{10}(L_2/L_{2,0})$',title=f'{titles[tag]}\n'+rf'$\rho={rho:.2f}$; {design} design')
                stem=f'{prefix}_{tag}';save_figure(fig,out,stem)
                atlas.append({'file':stem,'caption':f'{titles[tag]}. Except naive_A, normalized on U (detectable in at least one band at the deepest limit). Hatches mark regions without two resolved luminosities; these are masked rather than extrapolated. No smoothing. Rho {rho}; {design}.'})
            # Signed discrepancy only in the identifiable two-luminosity region.
            d=(qfit-qtrue).reshape(model.shape)[::-1,::-1][finite]/areas
            resolved=model.common_detected.reshape(model.shape)[::-1,::-1][finite]
            fig,ax=plt.subplots(figsize=(6.1,5.25))
            im=ax.pcolormesh(x1,x2,np.ma.masked_where(~resolved,d).T,norm=Normalize(vmin=-.25,vmax=.25),shading='flat',rasterized=True)
            fig.colorbar(im,ax=ax,pad=.025).set_label(r'Estimated minus true density [dex$^{-2}$]')
            ax.axvline(cutoff[0],ls=':',lw=1.1);ax.axhline(cutoff[1],ls=':',lw=1.1)
            ax.set(xlabel=r'$\log_{10}(L_1/L_{1,0})$',ylabel=r'$\log_{10}(L_2/L_{2,0})$',title=f'Reconstruction error in resolved cells\n'+rf'$\rho={rho:.2f}$; {design} design')
            stem=f'{prefix}_residual_resolved';save_figure(fig,out,stem)
            atlas.append({'file':stem,'caption':'Signed Q-density error, shown only where both luminosities are resolved by the deepest field. Masked cells have no displayed pointwise error. The fixed density range is shared across designs and correlations.'})
            # Marginals conditional on the common, fully resolved rectangle W1,W2 >= -0.5.
            low=np.r_[-np.inf,we];selected=low>=-.5;centers=(x1[:-1]+x1[1:])/2
            for b in (0,1):
                fig,ax=plt.subplots(figsize=(6.2,4.5))
                for label,m in [('True conditional LF',truth),('D-completion EM',fit),('No-D ablation',noD),('A-only histogram',naive)]:
                    mx=m.reshape(model.shape)[::-1,::-1].copy()
                    mask=selected[:,None]&selected[None,:];mx[~mask]=0;mx/=mx.sum()
                    marginal=mx.sum(axis=1-b)[1:-1]/np.diff((x1,x2)[b])
                    ax.stairs(marginal,(x1,x2)[b],label=label,linewidth=1.8 if label.startswith('True') else 1.15,linestyle='--' if label.startswith('True') else '-')
                ax.set(xlim=(config['mean_log10L'][b]-.5*config['sigma_dex'][b],(x1,x2)[b][-1]),xlabel=rf'$\log_{{10}}(L_{b+1}/L_{{{b+1},0}})$',ylabel=r'Conditional density [dex$^{-1}$]',title=f'Luminosity marginal in common resolved region\n'+rf'$\rho={rho:.2f}$; {design}')
                ax.legend(fontsize=8,frameon=False);ax.grid(alpha=.2)
                stem=f'{prefix}_marginal{b+1}';save_figure(fig,out,stem)
                atlas.append({'file':stem,'caption':'Marginal density conditional on both standardized log luminosities >= -0.5. This is not an extrapolated full-parent marginal.'})
            # Same-likelihood continuum: parent bright probabilities depend on unknown theta.
            a,qq=model.probabilities(qfit);theta_max=1-.5/min(qq)
            if theta_max < -1e-10:raise ValueError('Visible representative violates half-tail class')
            theta_max=max(0.,theta_max)
            ref_inv=fit[model.invisible].sum()
            R=np.zeros(model.H);R[model.invisible]=fit[model.invisible]/ref_inv
            thetas=[0.,theta_max/2,theta_max]
            public=json.loads((root/'public_counts'/design/f'{case}.json').read_text());counts=[np.asarray(n) for n in public['counts']]
            def ll(m):
                aa,qt=model.probabilities(m)
                return sum(float(ns[ns>0]@(np.log(ai[ns>0])-np.log(qj))) for ns,ai,qj in zip(counts,aa,qt))
            ll0=ll(qfit);ds=[];curves=[];families=[];maxlaw=0.
            low=np.r_[-np.inf,we];mask=np.isfinite(low)&(low>=-.5)&(low<=3)
            j0=int(np.flatnonzero(low==0)[0]);xplot=config['mean_log10L'][0]+config['sigma_dex'][0]*low[mask]
            fig,ax=plt.subplots(figsize=(6.3,4.7))
            tt=bright_surface(truth.reshape(model.shape)[::-1,::-1])[mask,j0]
            ax.plot(xplot,tt,ls='--',lw=2.,label='True parent probability')
            for theta in thetas:
                mm=(1-theta)*qfit+theta*R;families.append(mm)
                ai,qi=model.probabilities(mm)
                maxlaw=max(maxlaw,max(float(np.abs(aii/qii-ajj/qjj).max()) for aii,qii,ajj,qjj in zip(ai,qi,a,qq)))
                curve=bright_surface(mm.reshape(model.shape)[::-1,::-1])[mask,j0]
                ax.plot(xplot,curve,label=rf'Same record law: $\theta={theta:.3f}$')
                curves.append(curve);ds.append(ll(mm)-ll0)
            ax.set(xlabel=r'$x_1=\log_{10}(L_1/L_{1,0})$',ylabel=r'$P(X_1\geq x_1,\ X_2\geq 10.2)$',title=rf'Parent normalization remains undetermined: $\rho={rho:.2f}$'+'\nSame-likelihood family, NOT a confidence band')
            ax.legend(fontsize=8,frameon=False);ax.grid(alpha=.2)
            stem=f'{prefix}_same_likelihood';save_figure(fig,out,stem)
            atlas.append({'file':stem,'caption':'Three parent distributions sharing the fitted record law and satisfying the declared tail bounds. This is an observational-equivalence illustration, not a 95% confidence interval or a global-optimality certificate.'})
            fm={'thetas':thetas,'reference_theta':float(ref_inv),'theta_max':theta_max,'true_theta':float(truth[model.invisible].sum()),'max_record_probability_difference':maxlaw,'log_likelihood_differences':ds,
                'min_selection_probabilities':[float(min(model.probabilities(m)[1])) for m in families]}
            family_metadata[prefix]=fm
            write_csv(root/'curve_data'/f'{prefix}_same_likelihood.csv',[{'x1':float(xv),'true_parent':float(tv),**{f'family_{j}':float(curves[j][i]) for j in range(3)}} for i,(xv,tv) in enumerate(zip(xplot,tt))])
            np.savez_compressed(root/'curve_data'/f'{prefix}_family_mass.npz',masses=np.array(families),thetas=thetas)
            for metric,title,ylabel in [('bright_probability_RMSE_Q','Accuracy of observable bivariate probabilities','RMSE of joint bright probability'),('resolved_bin_L1_Q','Recovery of resolved luminosity-bin masses',r'$\sum_{h\in R}|\widehat Q_h-Q_{0,h}|$')]:
                fig,ax=plt.subplots(figsize=(6.1,4.7))
                for method in LABELS:
                    ns=[];means=[];lowerq=[];upperq=[]
                    for size in config['parent_sizes']:
                        vals=np.array([float(r[metric]) for r in rows if float(r['rho'])==rho and r['design']==design and int(r['parent_n'])==size and r['method']==method])
                        ns.append(size);means.append(vals.mean());lowerq.append(np.quantile(vals,.1));upperq.append(np.quantile(vals,.9))
                    ax.plot(ns,means,'o-',label=LABELS[method]);ax.fill_between(ns,lowerq,upperq,alpha=.12)
                ax.set(xscale='log',yscale='log',xlabel='Number of parent objects (simulation only)',ylabel=ylabel,title=f'{title}\n'+rf'$\rho={rho:.2f}$; {design}; 20 independent parent catalogues')
                ax.legend(fontsize=8,frameon=False);ax.grid(alpha=.2)
                stem=f'{prefix}_{metric}';save_figure(fig,out,stem)
                atlas.append({'file':stem,'caption':'Mean error over 20 realizations; shaded regions are empirical 10th--90th percentiles, NOT confidence bands. The same parent catalogue is reused across the two depth designs.'})
    # Actual public records, plotted as bin centers and upper-limit arrows.
    rho=config['rho'][0];case=f"rho{rho:.2f}_N{config['representative']['parent_size']}_r000"
    model=GridModel.from_dict(json.loads((root/'models/shallow.json').read_text()))
    with np.load(root/'latent_evaluation_only'/f'{case}.npz') as v:
        x=v['x'].copy();field=v['field'].copy()
    inc,delta,labels,idx=model.record_catalogue(-x,field)
    xe=[-e[::-1] for e in model.z_edges]
    def centers(b,indices):
        # The original values are not used to position detected symbols.
        zedges=np.r_[-np.inf,model.z_edges[b],np.inf]
        vals=[]
        for k in indices:
            if not np.isfinite(zedges[k:k+2]).all():vals.append(np.nan)
            else:vals.append(-.5*(zedges[k]+zedges[k+1]))
        return np.array(vals)
    fig,ax=plt.subplots(figsize=(6.6,5.3))
    for name,a,b,cap in [('A: both recorded',1,1,350),('B: band 2 upper limit',1,0,80),('C: band 1 upper limit',0,1,80)]:
        ids=np.flatnonzero((delta[:,0]==a)&(delta[:,1]==b))[:cap]
        xx=centers(0,labels[ids,2]) if a else -model.limits_z[field[ids],0]
        yy=centers(1,labels[ids,3]) if b else -model.limits_z[field[ids],1]
        ok=np.isfinite(xx)&np.isfinite(yy);xx,yy=xx[ok],yy[ok]
        if a and b:ax.scatter(xx,yy,s=9,alpha=.45,label=name)
        elif a:ax.errorbar(xx,yy,yerr=.18,uplims=True,fmt='.',markersize=3,alpha=.55,linestyle='none',label=name)
        else:ax.errorbar(xx,yy,xerr=.18,xuplims=True,fmt='.',markersize=3,alpha=.55,linestyle='none',label=name)
    ax.set(xlim=(xe[0][0],xe[0][-1]),ylim=(xe[1][0],xe[1][-1]),xlabel=r'$\log_{10}(L_1/L_{1,0})$',ylabel=r'$\log_{10}(L_2/L_{2,0})$',title='The catalogue supplied to the estimator\nBin records and upper limits; D objects absent')
    ax.legend(fontsize=8,frameon=False,loc='upper left');ax.grid(alpha=.15)
    save_figure(fig,out,'observed_catalogue')
    atlas.insert(0,{'file':'observed_catalogue','caption':'A fixed subset of public A/B/C records in the shallow Gaussian rho=0.55 catalogue. Detected markers use recorded-bin centers solely for plotting. Arrow lengths are graphical, not measured uncertainties. Hidden luminosities and D objects are not plotted.'})
    save_json(root/'same_likelihood_diagnostics.json',family_metadata)
    save_json(root/'figure_manifest.json',atlas)
    return atlas
