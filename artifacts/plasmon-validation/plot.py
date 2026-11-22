"""Plot measured public-solver CSVs; analytic lines are independent references."""
from pathlib import Path
import os,sys,glob,csv,json
os.environ.setdefault('MPLCONFIGDIR','/tmp/plasmon-validation/mpl-cache')
# This workstation has plotting packages in Nix store paths. In a conventional
# Python environment, installing numpy/matplotlib suffices and this is a no-op.
sys.path.extend(glob.glob('/nix/store/*/lib/python3.13/site-packages'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parent
with (ROOT/'roots.csv').open() as f: rows=list(csv.DictReader(f))
with (ROOT/'characteristics.csv').open() as f: maps=list(csv.DictReader(f))
with (ROOT/'map_roots.csv').open() as f: centers={r['case']:r for r in csv.DictReader(f)}
def group(case):
    rs=[r for r in rows if r['case']==case]
    return {k:np.array([float(r[k]) for r in rs]) for k in rs[0] if k not in ['case','classification']}
b=group('beta');w=group('omega');eta=.02
blue='#126f86';orange='#c66b23';gray='#737b85';purple='#9562ac'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white','axes.grid':True,'grid.alpha':.17,'legend.fontsize':9,'legend.frameon':False})
def eps(u):return 1-1/(u*(u+1j*eta))
def exact_b(u):return u*np.sqrt(eps(u)/(1+eps(u)))
def exact_w(v):
    u=np.sqrt(2*v*v/(1+2*v*v+np.sqrt(1+4*v**4))).astype(complex)-1j*eta/4
    for _ in range(25):
        u-=(u**4+1j*eta*u**3-(1+2*v*v)*u*u-2j*eta*v*v*u+v*v)/(4*u**3+3j*eta*u*u-2*(1+2*v*v)*u-2j*eta*v*v)
    return u
def save(fig,name):
    fig.savefig(ROOT/(name+'.png'),dpi=190,bbox_inches='tight')
    fig.savefig(ROOT/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(ROOT/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)
model=r'Vacuum | Drude metal: $\varepsilon_m=1-[u(u+0.02i)]^{-1}$, $u=\omega/\omega_p$, $k_p=10^5\,\mathrm{cm}^{-1}$'
u=np.linspace(.15,.695,750); ba=exact_b(u)
v=np.linspace(.2,3,750); wa=exact_w(v)
vl=np.linspace(.0,3.,750); ul=np.sqrt(2*vl*vl/(1+2*vl*vl+np.sqrt(1+4*vl**4)))
fig,ax=plt.subplots(2,2,figsize=(12,8.3),layout='constrained')
fig.suptitle('A lossy interface plasmon, solved in both complex cases',fontsize=18,fontweight='bold')
fig.supxlabel(model,fontsize=10)
a=ax[0,0];a.plot(ba.real,u,color=orange,lw=2,label='Independent TM relation');a.scatter(b['re'],b['parameter'],s=25,facecolors='white',edgecolors=orange,lw=1.3,zorder=4,label='lamina-mode accepted roots');a.plot(vl,ul,'--',c=gray,lw=1,label='Lossless reference');a.plot([0,.75],[0,.75],':',c=gray,lw=1,label='Vacuum light line');a.set(xlabel=r'$\mathrm{Re}(\beta)/k_p$',ylabel=r'$\omega/\omega_p$ (real)',title='A   Real frequency → complex wavevector',xlim=(0,3),ylim=(0,.76));a.legend(loc='lower right')
a=ax[0,1];a.plot(v,wa.real,color=blue,lw=2,label='Independent quartic root');a.scatter(w['parameter'],w['re'],s=25,facecolors='white',edgecolors=blue,lw=1.3,zorder=4,label='lamina-mode accepted roots');a.plot(vl,ul,'--',c=gray,lw=1,label='Lossless reference');a.axhline(1/np.sqrt(2),c=gray,ls=':',lw=1);a.text(1.2,.716,r'$\omega_{sp}/\omega_p=1/\sqrt{2}$',fontsize=9,color=gray);a.set(xlabel=r'$\beta/k_p$ (real)',ylabel=r'$\mathrm{Re}(\omega)/\omega_p$',title='B   Real wavevector → complex frequency',xlim=(0,3),ylim=(0,.76));a.legend(loc='lower right')
a=ax[1,0];a.plot(u,ba.imag,color=orange,lw=2);a.scatter(b['parameter'],b['im'],s=25,facecolors='white',edgecolors=orange,lw=1.3,zorder=4);a.set(xlabel=r'$\omega/\omega_p$ (real)',ylabel=r'$\mathrm{Im}(\beta)/k_p$',title='C   Positive spatial attenuation');a.text(.03,.93,r'$\exp(i\beta z)\;\Rightarrow\;\mathrm{Im}(\beta)>0$',transform=a.transAxes,va='top',color=orange)
a=ax[1,1];a.plot(v,-wa.imag,color=blue,lw=2);a.scatter(w['parameter'],-w['im'],s=25,facecolors='white',edgecolors=blue,lw=1.3,zorder=4);a.axhline(eta/2,c=gray,ls='--',lw=1,label=r'High-$\beta$ limit: $\gamma/2$');a.set(xlabel=r'$\beta/k_p$ (real)',ylabel=r'$-\mathrm{Im}(\omega)/\omega_p$',title='D   Positive temporal decay rate',ylim=(0,.011));a.text(.03,.9,r'$\exp(-i\omega t)\;\Rightarrow\;\mathrm{Im}(\omega)<0$',transform=a.transAxes,va='top',color=blue);a.legend(loc='lower right')
save(fig,'dispersion')
fig,axes=plt.subplots(1,2,figsize=(12.3,5.5),layout='constrained')
fig.suptitle('The solver characteristic actually vanishes',fontsize=18,fontweight='bold')
fig.supxlabel('Checked characteristic F from the public API; contours cross at the accepted numerical root.',fontsize=10)
for a,case,title,xlabel,ylabel in zip(axes,['beta','omega'],[r'Real $\omega=0.60\,\omega_p$',r'Real $\beta=1.00\,k_p$'],[r'$\mathrm{Re}(\beta)/k_p$',r'$\mathrm{Re}(\omega)/\omega_p$'],[r'$\mathrm{Im}(\beta)/k_p$',r'$\mathrm{Im}(\omega)/\omega_p$']):
    rs=[r for r in maps if r['case']==case]
    x=np.array([float(r['x']) for r in rs]).reshape(61,81);y=np.array([float(r['y']) for r in rs]).reshape(61,81)
    z=np.array([complex(float(r['real']),float(r['imag'])) for r in rs]).reshape(61,81)
    im=a.pcolormesh(x,y,np.log10(np.maximum(abs(z),1e-15)),cmap='magma',vmin=-4,vmax=-.25,shading='auto',rasterized=True)
    a.contour(x,y,z.real,levels=[0],colors=['#38dfdd'],linewidths=1.7)
    a.contour(x,y,z.imag,levels=[0],colors=['#9ade56'],linewidths=1.7)
    r=centers[case]
    a.scatter([float(r['re'])],[float(r['im'])],s=130,facecolors='none',edgecolors='white',lw=2,zorder=8)
    a.scatter([float(r['exact_re'])],[float(r['exact_im'])],marker='+',s=60,c='black',lw=1.5,zorder=9)
    a.set(xlabel=xlabel,ylabel=ylabel,title=title);a.grid(False)
    fig.colorbar(im,ax=a,label=r'$\log_{10}|F/k_p|$',shrink=.85)
handles=[Line2D([],[],color='#38bdbb',lw=2,label=r'Re $F=0$'),Line2D([],[],color='#7eaf37',lw=2,label=r'Im $F=0$'),Line2D([],[],color='black',marker='o',markerfacecolor='none',lw=0,label='Accepted numerical root'),Line2D([],[],color='black',marker='+',lw=0,label='Independent analytic root')]
fig.legend(handles=handles,loc='outside lower center',ncol=4,bbox_to_anchor=(.5,-.055))
save(fig,'characteristic_maps')
fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
fig.suptitle('Physicality checks beyond finding a small determinant',fontsize=18,fontweight='bold')
fig.supxlabel('Explicit decay-branch references; these checks do not constitute global mode coverage or QNM normalization.',fontsize=10)
for a,data,label,color in zip(axes[0],[b,w],[r'Real $\omega$; complex $\beta$',r'Real $\beta$; complex $\omega$'],[orange,blue]):
    x=data['parameter'];err=np.abs((data['re']+1j*data['im'])-(data['exact_re']+1j*data['exact_im']))/np.abs(data['exact_re']+1j*data['exact_im'])
    a.semilogy(x,np.maximum(err,1e-17),'o-',c=color,ms=3,label='Relative root error');a.semilogy(x,np.maximum(data['dispersion_error'],1e-17),'s--',c=gray,ms=3,label='TM boundary mismatch / kp');a.set(title=label,ylabel='Dimensionless error',xlabel=r'Prescribed real $\omega/\omega_p$' if data is b else r'Prescribed real $\beta/k_p$');a.legend()
for a,data,label,color in zip(axes[1],[b,w],[r'Real $\omega$: exterior decay',r'Real $\beta$: exterior decay'],[orange,blue]):
    x=data['parameter'];a.plot(x,data['kappa_left_im'],c=color,lw=2,label='Vacuum side');a.plot(x,data['kappa_right_im'],c=gray,lw=2,label='Metal side');a.set(title=label,ylabel=r'$\mathrm{Im}(\kappa)/k_p$  (must be positive)',xlabel=r'Prescribed real $\omega/\omega_p$' if data is b else r'Prescribed real $\beta/k_p$',ylim=(0,None));a.legend()
save(fig,'physicality_checks')
summary={}
for case,data in [('complex_beta',b),('complex_frequency',w)]:
    root=data['re']+1j*data['im'];exact=data['exact_re']+1j*data['exact_im']
    summary[case]={'accepted_points':len(root),'max_relative_root_error':float(np.max(abs(root-exact)/abs(exact))),'max_boundary_mismatch_over_kp':float(max(data['dispersion_error'])),'max_characteristic_residual_over_kp':float(max(data['residual'])),'max_scatter_transfer_value_difference_over_kp':float(max(data['backend_relative_error'])),'min_dielectric_decay_over_kp':float(min(data['kappa_left_im'])),'min_metal_decay_over_kp':float(min(data['kappa_right_im'])),'te_modes_at_sampled_checks':int(sum(data['te_accepted']))}
summary['model']={'epsilon_infinity':1,'epsilon_dielectric':1,'gamma_over_omega_p':eta,'kp_inverse_cm':1e5,'fitted_material':False,'te_checks_per_case':5,'frequency_reference':'explicit complex; built-in radiation qualification remains unspecified','branch_tracking':'independent local point searches; analytic solution locates valid benchmark charts'}
(ROOT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
html='''<!doctype html><html><meta charset="utf-8"><title>Planar plasmon validation</title><style>body{font:17px system-ui;max-width:1150px;margin:40px auto;padding:0 24px;color:#192933;line-height:1.55}img{width:100%}code,pre{background:#f1f4f5}pre{padding:18px;overflow:auto}h1,h2{line-height:1.2}</style><h1>Planar plasmon: numerical and physical validation</h1><p>Two half-spaces: vacuum | passive Drude metal. This is an illustrative analytical metal, not a fit to gold or silver. Time convention exp(−iωt); propagation exp(+iβz). The complex-wavevector and complex-frequency plots are different real-axis slices of the same dispersion surface, not interchangeable inverses at finite loss.</p><p><b>What was run:</b> 33 public APM discovery/refinement solves at real frequency, 33 at real wavevector, two additional solves for the characteristic maps, plus five local TE controls per case. Analytical roots place the small benchmark charts and select exterior branches; the markers are independent numerically accepted solver roots, not values copied from the formula. These charts do not demonstrate blind global discovery.</p><img src="dispersion.png"><img src="characteristic_maps.png"><img src="physicality_checks.png"><h2>Independent equations and acceptance checks</h2><p>Let u=ω/ωp and v=β/kp. εm=1−1/[u(u+iη)], η=0.02. At real u, v=u sqrt[εm/(1+εm)], with the passive forward sign. At real v the physical root obeys u⁴+iηu³−(1+2v²)u²−2iηv²u+v²=0, selected continuously from the positive lossless bound branch. The unsquared TM boundary condition κd/εd+κm/εm=0 is checked independently at each accepted point, with Im κd, Im κm positive. Transfer2 and Scatter2 checked values are compared at every root. The two half-spaces are used unchanged in both numerical paths.</p><h2>Measured results</h2><pre>'''+json.dumps(summary,indent=2)+'''</pre><h2>Interpretation limits</h2><p>The eigenvalues and exterior decay support the physical SPP interpretation for this benchmark. The built-in named dielectric classifier deliberately marks negative-permittivity exteriors unsupported, and explicit complex anchors retain unspecified built-in radiation provenance. We checked the bound physical branch independently; we did not manufacture a verified-positive-dielectric anchor, classify metal modes with the paper's dielectric labels, compute fields/QNM normalization, certify arbitrary materials or establish global mode coverage. TE controls only test the listed local windows.</p><h2>Reproduce</h2><p>From the lamina-core repository, with the Rust development environment loaded:</p><pre>cargo run --release --manifest-path artifacts/plasmon-validation/Cargo.toml --locked --offline -- artifacts/plasmon-validation\npython artifacts/plasmon-validation/plot.py</pre><p>Plotting requires NumPy and Matplotlib. Raw CSVs, Rust/Python source, the local Cargo manifest/lockfile and PNG/PDF/SVG figures are included. The implementation uses public lamina-mode API imports. No solver/library source was changed.</p><p>Benchmark background: <a href="https://doc.comsol.com/6.3/doc/com.comsol.help.models.woptics.metal_air_surface_plasmon_polariton/metal_air_surface_plasmon_polariton.html">COMSOL's metal–dielectric interface derivation</a> (TM boundary condition, evanescence and TE exclusion).</p></html>'''
(ROOT/'report.html').write_text(html)
