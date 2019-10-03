import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import warnings
warnings.filterwarnings('ignore')


T_edema = 10. #edema half-life
T_pal = 17. #palladium
T_iod = 60. #iodine
T_ces = 9.7 #cesium

#----------------------------------
#determine optimal post-treatment imaging time for three isotopes given edema half-life T_edema
if True:
	def get_T_bar(T_isotope, T_edema):
		l_edema = np.log(2) / T_edema
		l_isotope = np.log(2) / T_isotope
		T_bar = (1 / l_edema) * np.log((l_isotope + l_edema) / l_isotope) #using average size, robust to isotope uncertainty
		return T_bar
	print('')
	print('Cesium optimal treatment time: ' + str(get_T_bar(T_ces, T_edema)))
	print('Palladium optimal treatment time: ' + str(get_T_bar(T_pal, T_edema)))
	print('Iodine optimal treatment time: ' + str(get_T_bar(T_iod, T_edema)))
	print('')
#----------------------------------
#----------------------------------
#----------------------------------

#BED dose calculations according to eq 7b with adjustable parameters
#----------------------------------
if True:
	def BED(D, T_isotope=17., alpha=0.15, beta=0.05, Tp=42., T_repair=0.27): #T_repair in hours, else in days 
		mu = np.log(2) / (T_repair / 24)
		Teff = 1.44*T_isotope*np.log(alpha*D*Tp / T_isotope) #find effective time
		l = np.log(2) / T_isotope
		D_eff = (D)*(1 - np.exp(-l*Teff)) #compute total dose delivered by effective time
		D0_dot = l*D
		bracket_term = 1 - np.exp(-2*l*Teff) - (2*l*(1 - np.exp(-(mu + l)*Teff)) / (mu + l))
		other_term = (beta / alpha)*(D0_dot / (mu - l))*(1 / (1 - np.exp(-l*Teff)))
		RE = 1 + bracket_term*other_term
		return D_eff*RE - np.log(2*(Teff / alpha*Tp))


	fig, ax = plt.subplots()
	plt.subplots_adjust(left=0.25, bottom=0.25)
	D = np.linspace(40, 150, 300) #desired point doses from our isotope - > BED(D): 'delivered' point doses
	alpha0 = 0.15
	beta0 = 0.05
	Tp0 = 42.
	T_repair0 = 0.27 
	bed = BED(D, alpha=alpha0, beta = beta0, Tp=Tp0, T_repair = T_repair0)
	a, = plt.plot(D, BED(D, T_isotope=T_ces, alpha=alpha0, beta = beta0, Tp=Tp0, T_repair = T_repair0), 'blue', label='Cesium')
	b, = plt.plot(D, BED(D, T_isotope=T_pal, alpha=alpha0, beta = beta0, Tp=Tp0, T_repair = T_repair0), 'orange', label='Palladium')
	c, = plt.plot(D, BED(D, T_isotope=T_iod, alpha=alpha0, beta = beta0, Tp=Tp0, T_repair = T_repair0), 'red', label='Iodine')
	d, = plt.plot(D, D, 'k--', label='D = BED')
	ax.margins(x=0)

	axalpha = plt.axes([0.25, 0.02, 0.65, 0.01], facecolor='lightgoldenrodyellow')
	axbeta = plt.axes([0.25, 0.05, 0.65, 0.01], facecolor='lightgoldenrodyellow')
	axTp = plt.axes([0.25, 0.08, 0.65, 0.01], facecolor='lightgoldenrodyellow')
	axT_repair = plt.axes([0.25, 0.11, 0.65, 0.01], facecolor='lightgoldenrodyellow')

	###parameter ranges from Carlson et al confidence intervals [232, 233]
	salpha = Slider(axalpha, 'Alpha', 0.09, 0.35, valinit=alpha0, valstep=0.005)
	sbeta = Slider(axbeta, 'Beta', 0.05, 0.085, valinit=beta0, valstep=0.001)
	sTp = Slider(axTp, 'Tp [d]', 10., 60., valinit=Tp0, valstep=1.)
	sT_repair = Slider(axT_repair, 'T_repair [hr]', 0.25, 11., valinit=T_repair0, valstep=0.05) 

	def update(val):
		alpha = salpha.val
		beta = sbeta.val
		Tp = sTp.val
		T_repair = sT_repair.val
		a.set_ydata(BED(D, T_isotope=T_ces, alpha=alpha, beta=beta, Tp=Tp, T_repair=T_repair))
		b.set_ydata(BED(D, T_isotope=T_pal, alpha=alpha, beta=beta, Tp=Tp, T_repair=T_repair))
		c.set_ydata(BED(D, T_isotope=T_iod, alpha=alpha, beta=beta, Tp=Tp, T_repair=T_repair))
		fig.canvas.draw_idle()

	salpha.on_changed(update)
	sbeta.on_changed(update)
	sTp.on_changed(update)
	sT_repair.on_changed(update)

	ax.axes.set_title('Dale BED (Eq 7)')
	ax.axes.set_xlabel('Nominal Dose')
	ax.axes.set_ylabel('BED')
	ax.legend()	

	plt.show()

