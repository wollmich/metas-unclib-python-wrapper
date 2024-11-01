# Michael Wollensack METAS - 22.01.2019 - 22.10.2024

import os as _os
import sys as _sys
import numpy as np

_py2 = _sys.version_info[0] == 2
if _py2:
	_pm = " +/- "
else:
	_pm = u" \u00b1 "

_unclib_path = _os.path.dirname(_os.path.realpath(__file__))

import clr as _clr
_clr.AddReference("System.Windows.Forms")
_clr.AddReference(_os.path.join(_unclib_path, "Metas.UncLib.Core.dll"))
_clr.AddReference(_os.path.join(_unclib_path, "Metas.UncLib.DistProp.dll"))
_clr.AddReference(_os.path.join(_unclib_path, "Metas.UncLib.LinProp.dll"))
_clr.AddReference(_os.path.join(_unclib_path, "Metas.UncLib.MCProp.dll"))
_clr.AddReference(_os.path.join(_unclib_path, "Metas.UncLib.Optimization.dll"))

if _sys.platform == 'win32':
	from ctypes import windll as _windll
	_windll.shcore.SetProcessDpiAwareness(1)

from System import Array as _Array
from System import Byte as _Byte
from System.Windows.Forms import Form as _Form
from System.Windows.Forms import DockStyle as _DockStyle
from System.Drawing import Size as _Size
from Metas.UncLib.Core import Number as _Number
from Metas.UncLib.Core import Real as _Real
from Metas.UncLib.Core import Complex as _Complex
from Metas.UncLib.Core import Const2014 as _Const2014
from Metas.UncLib.Core import Const2014_90 as _Const2014_90
from Metas.UncLib.Core import Const2018 as _Const2018
from Metas.UncLib.Core import Math as _Math
from Metas.UncLib.Core import SplineBoundary
from Metas.UncLib.Core.Ndims import RealNArray as _RealNArray
from Metas.UncLib.Core.Ndims import ComplexNArray as _ComplexNArray
from Metas.UncLib.Core.Ndims import RealLinAlg as _RealLinAlg
from Metas.UncLib.Core.Ndims import ComplexLinAlg as _ComplexLinAlg
from Metas.UncLib.Core.Ndims import RealNumLib as _RealNumLib
from Metas.UncLib.Core.Ndims import ComplexNumLib as _ComplexNumLib
from Metas.UncLib.Core.Unc import InputId as _InputId
from Metas.UncLib.Core.Unc import GenericUnc as _GenericUnc
from Metas.UncLib.Core.Unc import StandardNormalDistribution
from Metas.UncLib.Core.Unc import NormalDistribution
from Metas.UncLib.Core.Unc import StandardUniformDistribution
from Metas.UncLib.Core.Unc import UniformDistribution
from Metas.UncLib.Core.Unc import CurvilinearTrapezoidDistribution
from Metas.UncLib.Core.Unc import TrapezoidalDistribution
from Metas.UncLib.Core.Unc import TriangularDistribution
from Metas.UncLib.Core.Unc import ArcSineDistribution
from Metas.UncLib.Core.Unc import ExponentialDistribution
from Metas.UncLib.Core.Unc import GammaDistribution
from Metas.UncLib.Core.Unc import ChiSquaredDistribution
from Metas.UncLib.Core.Unc import StudentTDistribution
from Metas.UncLib.Core.Unc import StudentTFromSamplesDistribution
from Metas.UncLib.Core.Unc import RandomChoicesFromSamples
from Metas.UncLib.DistProp import UncNumber as _DistPropUncNumber
from Metas.UncLib.DistProp import UncList as _DistPropUncList
from Metas.UncLib.DistProp.Misc import Global as _DistPropGlobal
from Metas.UncLib.LinProp import UncNumber as _LinPropUncNumber
from Metas.UncLib.LinProp import UncList as _LinPropUncList
from Metas.UncLib.LinProp.Misc import Global as _LinPropGlobal
from Metas.UncLib.LinProp.Misc import DofModeType
from Metas.UncLib.LinProp.Misc import FromSamplesModeType
from Metas.UncLib.LinProp.Ndims import RealUncLinAlg as _RealUncLinAlg
from Metas.UncLib.LinProp.Ndims import ComplexUncLinAlg as _ComplexUncLinAlg
from Metas.UncLib.LinProp.Ndims import RealUncNumLib as _RealUncNumLib
from Metas.UncLib.LinProp.Ndims import ComplexUncNumLib as _ComplexUncNumLib
from Metas.UncLib.LinProp import NumericalFunctionDelegate as _NumericalFunctionDelegate
from Metas.UncLib.LinProp import UncNumerical as _UncNumerical
from Metas.UncLib.LinProp import UnknownIdDecoderDelegate as _UnknownIdDecoderDelegate
from Metas.UncLib.LinProp.Gui import GuiUncBudget as _GuiUncBudget
from Metas.UncLib.LinProp.Gui import GuiUncListBudget as _GuiUncListBudget
from Metas.UncLib.LinProp.Gui import HighDPIHelper as _HighDPIHelper
from Metas.UncLib.MCProp import UncNumber as _MCPropUncNumber
from Metas.UncLib.MCProp import UncList as _MCPropUncList
from Metas.UncLib.MCProp.Misc import Global as _MCPropGlobal
from Metas.UncLib.Optimization import Optimizer as _Optimizer
from Metas.UncLib.Optimization import ObjectiveFunction as _ObjectiveFunction
from Metas.UncLib.Optimization import Algorithm as OptimizerAlgorithm

_UncNumber = _LinPropUncNumber
_UncList = _LinPropUncList
_UncHelper = _GenericUnc[_UncList, _UncNumber]()

def use_distprop(maxlevel=1):
	global _UncNumber
	global _UncList
	global _UncHelper
	_UncNumber = _DistPropUncNumber
	_UncList = _DistPropUncList
	_UncHelper = _GenericUnc[_UncList, _UncNumber]()
	_DistPropGlobal.MaxLevel = int(maxlevel)

def use_linprop(dofmode=DofModeType.WelchSatterthwaite, fromsamplesmode=FromSamplesModeType.ExpandInputCovariance):
	global _UncNumber
	global _UncList
	global _UncHelper
	_UncNumber = _LinPropUncNumber
	_UncList = _LinPropUncList
	_UncHelper = _GenericUnc[_UncList, _UncNumber]()
	_LinPropGlobal.DofMode = dofmode
	_LinPropGlobal.FromSamplesMode = fromsamplesmode

def use_mcprop(n=100000):
	global _UncNumber
	global _UncList
	global _UncHelper
	_UncNumber = _MCPropUncNumber
	_UncList = _MCPropUncList
	_MCPropGlobal.n = int(n)
	_UncHelper = _GenericUnc[_UncList, _UncNumber]()

def _asnetobject(a, forcecomplex=False, forcecomplexifnegative=False):
	if isinstance(a, ufloat):
		if forcecomplex:
			b = ucomplex(a).net_object
		elif forcecomplexifnegative and get_fcn_value(a) < 0:
			b = ucomplex(a).net_object
		else:
			b = a.net_object
	elif isinstance(a, ucomplex):
		b = a.net_object
	elif isinstance(a, float):
		if forcecomplex:
			b = ucomplex(a).net_object
		elif forcecomplexifnegative and a < 0:
			b = ucomplex(a).net_object
		else:
			b = ufloat(a).net_object
	elif isinstance(a, complex):
		b = ucomplex(a).net_object
	else:
		b = _asnetnarray(a, forcecomplex, forcecomplexifnegative)
	return b

def _fromnetobject(a):
	if type(a) is _Complex[_UncNumber]:
		b = ucomplex(a)
	elif type(a) is _UncNumber:
		b = ufloat(a)
	elif type(a) is _Complex[_Number]:
		b = complex(float(a.Real().Value), float(a.Imag().Value))
	elif type(a) is _Number:
		b = float(a.Value)
	else:
		b = _fromnetnarray(a)
	return b

def _asunclist(a):
	b = _asnetnarray(a)
	c = _UncList()
	return c.op_Implicit(b)

def _asnetnarray(a, complexarray=False, complexarrayifnegative=False):
	b = np.asarray(a)
	s = b.shape
	s2 = _Array[int](s)
	if complexarray or iscomplexarray(b):
		c = [ucomplex(bi).net_object for bi in b.flatten('F')]
		d = _ComplexNArray[_UncNumber]()
	elif complexarrayifnegative and any([get_fcn_value(ufloat(bi)) < 0 for bi in b.flatten('F')]):
		c = [ucomplex(bi).net_object for bi in b.flatten('F')]
		d = _ComplexNArray[_UncNumber]()
	else:
		c = [ufloat(bi).net_object for bi in b.flatten('F')]
		d = _RealNArray[_UncNumber]()
	d.Init1dData(c, False)
	d.Reshape(s2)
	return d

def _fromnetnarray(a):
	s = a.size
	s2 = [si for si in s]
	d = a.data
	if type(a) is _ComplexNArray[_UncNumber]:
		d2 = [ucomplex(di) for di in d]
	elif type(a) is _RealNArray[_UncNumber]:
		d2 = [ufloat(di) for di in d]
	elif type(a) is _ComplexNArray[_Number]:
		d2 = [complex(float(di.Real().Value), float(di.Imag().Value)) for di in d]
	else:
		d2 = [float(di.Value) for di in d]
	d3 = np.asarray(d2)
	d4 = d3.reshape(s2, order='F')
	return d4

def _asnetnumbernarray(a):
	b = np.asarray(a)
	s = b.shape
	s2 = _Array[int](s)
	c = [_Number(float(bi)) for bi in b.flatten('F')]
	d = _RealNArray[_Number]()
	d.Init1dData(c, False)
	d.Reshape(s2)
	return d

def _asnetcomplexnumbernarray(a):
	b = np.asarray(a)
	s = b.shape
	s2 = _Array[int](s)
	c = [_Complex[_Number](_Number(complex(bi).real), _Number(complex(bi).imag)) for bi in b.flatten('F')]
	d = _ComplexNArray[_Number]()
	d.Init1dData(c, False)
	d.Reshape(s2)
	return d

def _fromnetjaggednumbermatrix(a):
	b = _RealNArray[_Number]()
	b.Init2dData(a)
	s = b.size
	s2 = [si for si in s]
	d = b.data
	d2 = [float(di.Value) for di in d]
	d3 = np.asarray(d2)
	d4 = d3.reshape(s2, order='F')
	return d4

def _asnetdoublearray(a):
	b = np.asarray(a)
	c = _Array[float]([float(bi) for bi in b.flatten('F')])
	return c

def _input_id_desc(id=None, desc=None):
	if id is None:
		id2 = _InputId.NewInputId()
	else:
		id2 = _InputId(id)
	if type(desc) is str:
		desc2 = desc
	else:
		desc2 = ''
	return id2, desc2


def unc_budget(a, format=None, name=None, infos=None):
	b = _asunclist(a)
	if len(b.data) == 1:
		c = _GuiUncBudget()
	else:
		c = _GuiUncListBudget()
	if format is not None:
		c.Format = format
	if name is None:
		name2 = ''
	else:
		name2 = name
	if infos is not None:
		c.Infos = [info + '\n' for info in infos]
	if len(b.data) == 1:
		c.Value = b.data[0]
	else:
		c.Values = b
	c.Dock = _DockStyle.Fill
	f = _Form()
	f.Text = name2
	f.Size = _Size(_HighDPIHelper.ScaleValueX(640), _HighDPIHelper.ScaleValueY(480))
	f.Controls.Add(c)
	f.ShowDialog()


def ufloatsystem(value, sys_inputs, sys_sensitivities):
	if _UncNumber is not _LinPropUncNumber:
		raise Exception("ufloatsystem supports only LinProp uncertainty objects")
	d = _UncNumber(0)
	value2 = float(value)
	sys_inputs2 = _asunclist(sys_inputs).data
	sys_sensitivities2 = _asnetdoublearray(sys_sensitivities)
	d.Init(value2, sys_inputs2, sys_sensitivities2)
	return ufloat(d)

def ufloatarray(values, covariance, idof=0.0, id=None, desc=None):
	v = _asnetnumbernarray(values)
	cv = _asnetnumbernarray(covariance)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.RealUncNArray(v, cv.Matrix, float(idof), id2, desc2)
	return _fromnetnarray(d)

def ucomplexarray(values, covariance, idof=0.0, id=None, desc=None):
	v = _asnetcomplexnumbernarray(values)
	cv = _asnetnumbernarray(covariance)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.ComplexUncNArray(v, cv.Matrix, float(idof), id2, desc2)
	return _fromnetnarray(d)

def ufloatfromsamples(samples, id=None, desc=None, p=0.95):
	s = _asnetnumbernarray(samples)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.RealUncNumberFromSamples(s.Vector, id2, desc2, p)
	return _fromnetobject(d)

def ucomplexfromsamples(samples, id=None, desc=None, p=0.95):
	s = _asnetcomplexnumbernarray(samples)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.ComplexUncNumberFromSamples(s.Vector, id2, desc2, p)
	return _fromnetobject(d)

def ufloatarrayfromsamples(samples, id=None, desc=None, p=0.95):
	s = _asnetnumbernarray(samples)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.RealUncNArrayFromSamples(s.Matrix, id2, desc2, p)
	return _fromnetnarray(d)

def ucomplexarrayfromsamples(samples, id=None, desc=None, p=0.95):
	s = _asnetcomplexnumbernarray(samples)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.ComplexUncNArrayFromSamples(s.Matrix, id2, desc2, p)
	return _fromnetnarray(d)

def ufloatfromrandomchoices(samples, id=None, desc=None):
	s = _asnetnumbernarray(samples)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.RealUncNumberFromRandomChoices(s.Vector, id2, desc2)
	return _fromnetobject(d)

def ucomplexfromrandomchoices(samples, id=None, desc=None):
	s = _asnetcomplexnumbernarray(samples)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.ComplexUncNumberFromRandomChoices(s.Vector, id2, desc2)
	return _fromnetobject(d)

def ufloatarrayfromrandomchoices(samples, id=None, desc=None):
	s = _asnetnumbernarray(samples)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.RealUncNArrayFromRandomChoices(s.Matrix, id2, desc2)
	return _fromnetnarray(d)

def ucomplexarrayfromrandomchoices(samples, id=None, desc=None):
	s = _asnetcomplexnumbernarray(samples)
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.ComplexUncNArrayFromRandomChoices(s.Matrix, id2, desc2)
	return _fromnetnarray(d)

def ufloatfromdistribution(distribution, id=None, desc=None):
	id2, desc2 = _input_id_desc(id, desc)
	d = _UncHelper.RealUncNumberFromDistribution(distribution, id2, desc2)
	return _fromnetobject(d)

def get_value(a):
	return _fromnetobject(_UncHelper.GetValue(_asnetobject(a)))

def get_fcn_value(a):
	return _fromnetobject(_UncHelper.GetFcnValue(_asnetobject(a)))

def get_stdunc(a):
	return _fromnetobject(_UncHelper.GetStdUnc(_asnetobject(a)))

def get_idof(a):
	return _fromnetobject(_UncHelper.GetIDof(_asnetobject(a)))

def get_moment(a, n):
	return _fromnetobject(_UncHelper.GetMoment(_asnetobject(a), int(n)))

def get_coverage_interval(a, p=0.95):
	a2 = _asunclist(a)
	r1 = _UncHelper.GetCoverageInterval(a2, float(p))
	r2 = _fromnetjaggednumbermatrix(r1)
	return r2

def get_covariance(a):
	a2 = _asunclist(a)
	r1 = _UncHelper.GetCovariance(a2)
	r2 = _fromnetjaggednumbermatrix(r1)
	return r2

def get_correlation(a):
	a2 = _asunclist(a)
	r1 = _UncHelper.GetCorrelation(a2)
	r2 = _fromnetjaggednumbermatrix(r1)
	return r2

def get_jacobi(a):
	a2 = _asunclist(a)
	r1 = _UncHelper.GetJacobi(a2)
	r2 = _fromnetjaggednumbermatrix(r1)
	return r2

def get_jacobi2(a, b):
	a2 = _asunclist(a)
	b2 = _asunclist(b)
	r1 = _UncHelper.GetJacobi2(a2, b2)
	r2 = _fromnetjaggednumbermatrix(r1)
	return r2

def get_unc_component(a, b):
	a2 = _asunclist(a)
	b2 = _asunclist(b)
	r1 = _UncHelper.GetUncComponent(a2, b2)
	r2 = _fromnetjaggednumbermatrix(r1)
	return r2

def iscomplex(d):
	return isinstance(d, (complex, _Complex[_UncNumber], ucomplex))

def iscomplexarray(a):
	return any([iscomplex(ai) for ai in np.asarray(a).flat])

def iszero(a):
	b = np.asarray(a)
	c = _asnetnarray(b)
	d = all([ci.IsZero() for ci in c.data])
	return d

class umath(object):
	@staticmethod
	def abs(a):
		return _fromnetobject(_asnetobject(a).Abs())

	@staticmethod
	def exp(a):
		return _fromnetobject(_asnetobject(a).Exp())

	@staticmethod
	def log(a):
		return _fromnetobject(_asnetobject(a, forcecomplexifnegative=True).Log())

	@staticmethod
	def log10(a):
		return _fromnetobject(_asnetobject(a, forcecomplexifnegative=True).Log10())

	@staticmethod
	def pow(a, b):
		return a**b

	@staticmethod
	def sqrt(a):
		return _fromnetobject(_asnetobject(a, forcecomplexifnegative=True).Sqrt())

	@staticmethod
	def sin(a):
		return _fromnetobject(_asnetobject(a).Sin())

	@staticmethod
	def cos(a):
		return _fromnetobject(_asnetobject(a).Cos())

	@staticmethod
	def tan(a):
		return _fromnetobject(_asnetobject(a).Tan())

	@staticmethod
	def asin(a):
		return _fromnetobject(_asnetobject(a).Asin())

	@staticmethod
	def acos(a):
		return _fromnetobject(_asnetobject(a).Acos())

	@staticmethod
	def atan(a):
		return _fromnetobject(_asnetobject(a).Atan())

	@staticmethod
	def arcsin(a):
		return _fromnetobject(_asnetobject(a).Asin())

	@staticmethod
	def arccos(a):
		return _fromnetobject(_asnetobject(a).Acos())

	@staticmethod
	def arctan(a):
		return _fromnetobject(_asnetobject(a).Atan())

	@staticmethod
	def sinh(a):
		return _fromnetobject(_asnetobject(a).Sinh())

	@staticmethod
	def cosh(a):
		return _fromnetobject(_asnetobject(a).Cosh())

	@staticmethod
	def tanh(a):
		return _fromnetobject(_asnetobject(a).Tanh())

	@staticmethod
	def asinh(a):
		return _fromnetobject(_asnetobject(a).Asinh())

	@staticmethod
	def acosh(a):
		return _fromnetobject(_asnetobject(a).Acosh())

	@staticmethod
	def atanh(a):
		return _fromnetobject(_asnetobject(a).Atanh())

	@staticmethod
	def arcsinh(a):
		return _fromnetobject(_asnetobject(a).Asinh())

	@staticmethod
	def arccosh(a):
		return _fromnetobject(_asnetobject(a).Acosh())

	@staticmethod
	def arctanh(a):
		return _fromnetobject(_asnetobject(a).Atanh())

	@staticmethod
	def real(a):
		return _fromnetobject(_asnetobject(a, True).Real())
	
	@staticmethod
	def imag(a):
		return _fromnetobject(_asnetobject(a, True).Imag())

	@staticmethod
	def conj(a):
		return _fromnetobject(_asnetobject(a, True).Conj())

	@staticmethod
	def conjugate(a):
		return _fromnetobject(_asnetobject(a, True).Conj())

	@staticmethod
	def angle(a, deg=False):
		b = _fromnetobject(_asnetobject(a, True).Angle())
		if deg:
			b = b * (180.0 / np.pi)
		return b

	@staticmethod
	def phase(a, deg=False):
		return umath.angle(a, deg)

	@staticmethod
	def ellipk(a):
		if (iscomplex(a) | iscomplexarray(a)):
			raise Exception("Input must be real")
		return _fromnetobject(_asnetobject(a).Ellipk())

	@staticmethod
	def ellipe(a):
		if (iscomplex(a) | iscomplexarray(a)):
			raise Exception("Input must be real")
		return _fromnetobject(_asnetobject(a).Ellipe())


class ulinalg(object):
	@staticmethod
	def dot(a, b):
		complexarray = iscomplexarray(a) or iscomplexarray(b)
		a2 = _asnetnarray(a, complexarray)
		b2 = _asnetnarray(b, complexarray)
		if complexarray:
			c2 = _ComplexLinAlg[_UncNumber].Dot(a2, b2)
		else:
			c2 = _RealLinAlg[_UncNumber].Dot(a2, b2)
		c = _fromnetnarray(c2)
		return c

	@staticmethod
	def inv(a):
		complexarray = iscomplexarray(a)
		a2 = _asnetnarray(a, complexarray)
		if complexarray:
			b2 = _ComplexLinAlg[_UncNumber].Inv(a2)
		else:
			b2 = _RealLinAlg[_UncNumber].Inv(a2)
		b = _fromnetnarray(b2)
		return b

	@staticmethod
	def det(a):
		complexarray = iscomplexarray(a)
		a2 = _asnetnarray(a, complexarray)
		if complexarray:
			b2 = _ComplexLinAlg[_UncNumber].Det(a2)
			b = ucomplex(b2)
		else:
			b2 = _RealLinAlg[_UncNumber].Det(a2)
			b = ufloat(b2)
		return b

	@staticmethod
	def solve(a, y):
		complexarray = iscomplexarray(a) or iscomplexarray(y)
		a2 = _asnetnarray(a, complexarray)
		y2 = _asnetnarray(y, complexarray)
		if complexarray:
			x2 = _ComplexLinAlg[_UncNumber].Solve(a2, y2)
		else:
			x2 = _RealLinAlg[_UncNumber].Solve(a2, y2)
		x = _fromnetnarray(x2)
		return x

	@staticmethod
	def lstsqrsolve(a, y):
		complexarray = iscomplexarray(a) or iscomplexarray(y)
		a2 = _asnetnarray(a, complexarray)
		y2 = _asnetnarray(y, complexarray)
		if complexarray:
			x2 = _ComplexLinAlg[_UncNumber].LstSqrSolve(a2, y2)
		else:
			x2 = _RealLinAlg[_UncNumber].LstSqrSolve(a2, y2)
		x = _fromnetnarray(x2)
		return x

	@staticmethod
	def weightedlstsqrsolve(a, y, w):
		complexarray = iscomplexarray(a) or iscomplexarray(y)
		a2 = _asnetnarray(a, complexarray)
		y2 = _asnetnarray(y, complexarray)
		w2 = _asnetnarray(w, False)
		if complexarray:
			x2 = _ComplexLinAlg[_UncNumber].WeightedLstSqrSolve(a2, y2, w2)
		else:
			x2 = _RealLinAlg[_UncNumber].WeightedLstSqrSolve(a2, y2, w2)
		x = _fromnetnarray(x2)
		return x

	@staticmethod
	def generallstsqrsolve(a, y, v):
		complexarray = iscomplexarray(a) or iscomplexarray(y)
		a2 = _asnetnarray(a, complexarray)
		y2 = _asnetnarray(y, complexarray)
		v2 = _asnetnarray(v, False)
		if complexarray:
			x2 = _ComplexLinAlg[_UncNumber].GeneralLstSqrSolve(a2, y2, v2)
		else:
			x2 = _RealLinAlg[_UncNumber].GeneralLstSqrSolve(a2, y2, v2)
		x = _fromnetnarray(x2)
		return x

	@staticmethod
	def lu(a):
		complexarray = iscomplexarray(a)
		a2 = _asnetnarray(a, complexarray)
		if complexarray:
			res = _ComplexLinAlg[_UncNumber].Lu(a2)
		else:
			res = _RealLinAlg[_UncNumber].Lu(a2)
		l = _fromnetnarray(res.l)
		u = _fromnetnarray(res.u)
		p = _fromnetnarray(res.p)
		return l, u, p

	@staticmethod
	def cholesky(a):
		complexarray = iscomplexarray(a)
		a2 = _asnetnarray(a, complexarray)
		if complexarray:
			l2 = _ComplexLinAlg[_UncNumber].Cholesky(a2)
		else:
			l2 = _RealLinAlg[_UncNumber].Cholesky(a2)
		l = _fromnetnarray(l2)
		return l

	@staticmethod
	def qr(a):
		complexarray = iscomplexarray(a)
		a2 = _asnetnarray(a, complexarray)
		if complexarray:
			res = _ComplexLinAlg[_UncNumber].Qr(a2)
		else:
			res = _RealLinAlg[_UncNumber].Qr(a2)
		q = _fromnetnarray(res[0])
		r = _fromnetnarray(res[1])
		return q, r

	@staticmethod
	def svd(a):
		complexarray = iscomplexarray(a)
		a2 = _asnetnarray(a, complexarray)
		if complexarray:
			res = _ComplexLinAlg[_UncNumber].Svd(a2)
		else:
			res = _RealLinAlg[_UncNumber].Svd(a2)
		u = _fromnetnarray(res[0])
		s = _fromnetnarray(res[1])
		v = _fromnetnarray(res[2])
		return u, s, v

	@staticmethod
	def eig(a, *args):
		if _UncNumber is not _LinPropUncNumber:
			raise Exception("Eig supports only LinProp uncertainty objects")
		b = [np.asarray(arg) for arg in args]
		b.insert(0, np.asarray(a))
		n1 = b[0].shape[0]
		n2 = b[0].shape[1]
		complexarray = any([iscomplexarray(bi) for bi in b])
		if complexarray or n1 != n2 or len(b) > 1:
			symmetric = False
		else:
			symmetric = all([iszero(bi - bi.T) for bi in b])
		complexeig = complexarray or not symmetric
		c = [_asnetnarray(bi, complexeig) for bi in b]
		if n1 != n2 and len(c) == 1:
			raise Exception("Matrix must be square")
		if len(c) == 1:
			if complexeig:
				res = _ComplexUncLinAlg[_UncNumber].NonsymmetricEig(c[0])
			else:
				res = _RealUncLinAlg[_UncNumber].SymmetricEig(c[0])
		else:
			res = _ComplexUncLinAlg[_UncNumber].NonLinearEig(c)
		v2 = _fromnetnarray(res[0])
		d2 = _fromnetnarray(res[1])
		v3 = v2[0:n2, :]
		d3 = d2
		if not complexarray and not symmetric and iszero(umath.imag(v3)) and iszero(umath.imag(d3)):
			v4 = umath.real(v3)
			d4 = umath.real(d3)
		else:
			v4 = v3
			d4 = d3
		return v4, d4


class unumlib(object):
	@staticmethod
	def polyfit(x, y, n):
		complexarray = iscomplexarray(x) or iscomplexarray(y)
		x2 = _asnetnarray(x, complexarray)
		y2 = _asnetnarray(y, complexarray)
		n2 = int(n)
		if complexarray:
			p2 = _ComplexNumLib[_UncNumber].PolyFit(x2, y2, n2)
		else:
			p2 = _RealNumLib[_UncNumber].PolyFit(x2, y2, n2)
		p = _fromnetnarray(p2)
		return p

	@staticmethod
	def polyval(p, x):
		complexarray = iscomplexarray(p) or iscomplexarray(x)
		p2 = _asnetnarray(p, complexarray)
		x2 = _asnetnarray(x, complexarray)
		if complexarray:
			y2 = _ComplexNumLib[_UncNumber].PolyVal(p2, x2)
		else:
			y2 = _RealNumLib[_UncNumber].PolyVal(p2, x2)
		y = _fromnetnarray(y2)
		return y

	@staticmethod
	def interpolation(x, y, n, xx):
		complexarray = iscomplexarray(y)
		x2 = _asnetdoublearray(x)
		y2 = _asnetnarray(y, complexarray)
		xx2 = _asnetdoublearray(xx)
		n2 = int(n)
		if complexarray:
			yy2 = _ComplexNumLib[_UncNumber].Interpolation(x2, y2, n2, xx2)
		else:
			yy2 = _RealNumLib[_UncNumber].Interpolation(x2, y2, n2, xx2)
		yy = _fromnetnarray(yy2)
		return yy

	@staticmethod
	def interpolation2(x, y, n, xx):
		complexarray = iscomplexarray(y)
		x2 = _asnetdoublearray(x)
		y2 = _asnetnarray(y, complexarray)
		xx2 = _asnetdoublearray(xx)
		n2 = int(n)
		if complexarray:
			yy2 = _ComplexNumLib[_UncNumber].Interpolation2(x2, y2, n2, xx2)
		else:
			yy2 = _RealNumLib[_UncNumber].Interpolation2(x2, y2, n2, xx2)
		yy = _fromnetnarray(yy2)
		return yy

	@staticmethod
	def splineinterpolation(x, y, xx, 
			startboundary=SplineBoundary.Natural_Spline, startderivativevalue=0.0,
			endboundary=SplineBoundary.Natural_Spline, endderivativevalue=0.0):
		complexarray = iscomplexarray(y)
		x2 = _asnetdoublearray(x)
		y2 = _asnetnarray(y, complexarray)
		xx2 = _asnetdoublearray(xx)
		if complexarray:
			sd = ucomplex(startderivativevalue).net_object
			ed = ucomplex(endderivativevalue).net_object
			yy2 = _ComplexNumLib[_UncNumber].SplineInterpolation(x2, y2, xx2, startboundary, sd, endboundary, ed)
		else:
			sd = ufloat(startderivativevalue).net_object
			ed = ufloat(endderivativevalue).net_object
			yy2 = _RealNumLib[_UncNumber].SplineInterpolation(x2, y2, xx2, startboundary, sd, endboundary, ed)
		yy = _fromnetnarray(yy2)
		return yy

	@staticmethod
	def splineinterpolation2(x, y, xx, 
			startboundary=SplineBoundary.Natural_Spline, startderivativevalue=0.0,
			endboundary=SplineBoundary.Natural_Spline, endderivativevalue=0.0):
		complexarray = iscomplexarray(y)
		x2 = _asnetdoublearray(x)
		y2 = _asnetnarray(y, complexarray)
		xx2 = _asnetdoublearray(xx)
		if complexarray:
			sd = ucomplex(startderivativevalue).net_object
			ed = ucomplex(endderivativevalue).net_object
			yy2 = _ComplexNumLib[_UncNumber].SplineInterpolation2(x2, y2, xx2, startboundary, sd, endboundary, ed)
		else:
			sd = ufloat(startderivativevalue).net_object
			ed = ufloat(endderivativevalue).net_object
			yy2 = _RealNumLib[_UncNumber].SplineInterpolation2(x2, y2, xx2, startboundary, sd, endboundary, ed)
		yy = _fromnetnarray(yy2)
		return yy

	@staticmethod
	def integrate(x, y, n):
		complexarray = iscomplexarray(y)
		x2 = _asnetdoublearray(x)
		y2 = _asnetnarray(y, complexarray)
		n2 = int(n)
		if complexarray:
			z2 = _ComplexNumLib[_UncNumber].Integrate(x2, y2, n2)
		else:
			z2 = _RealNumLib[_UncNumber].Integrate(x2, y2, n2)
		z = _fromnetnarray(z2)
		return z

	@staticmethod
	def integrate2(x, y, n):
		complexarray = iscomplexarray(y)
		x2 = _asnetdoublearray(x)
		y2 = _asnetnarray(y, complexarray)
		n2 = int(n)
		if complexarray:
			z2 = _ComplexNumLib[_UncNumber].Integrate2(x2, y2, n2)
			z = ucomplex(z2)
		else:
			z2 = _RealNumLib[_UncNumber].Integrate2(x2, y2, n2)
			z = ufloat(z2)
		return z

	@staticmethod
	def splineintegrate(x, y, 
			startboundary=SplineBoundary.Natural_Spline, startderivativevalue=0.0,
			endboundary=SplineBoundary.Natural_Spline, endderivativevalue=0.0):
		complexarray = iscomplexarray(y)
		x2 = _asnetdoublearray(x)
		y2 = _asnetnarray(y, complexarray)
		if complexarray:
			sd = ucomplex(startderivativevalue).net_object
			ed = ucomplex(endderivativevalue).net_object
			z2 = _ComplexNumLib[_UncNumber].SplineIntegrate(x2, y2, startboundary, sd, endboundary, ed)
		else:
			sd = ufloat(startderivativevalue).net_object
			ed = ufloat(endderivativevalue).net_object
			z2 = _RealNumLib[_UncNumber].SplineIntegrate(x2, y2, startboundary, sd, endboundary, ed)
		z = _fromnetnarray(z2)
		return z

	@staticmethod
	def splineintegrate2(x, y, 
			startboundary=SplineBoundary.Natural_Spline, startderivativevalue=0.0,
			endboundary=SplineBoundary.Natural_Spline, endderivativevalue=0.0):
		complexarray = iscomplexarray(y)
		x2 = _asnetdoublearray(x)
		y2 = _asnetnarray(y, complexarray)
		if complexarray:
			sd = ucomplex(startderivativevalue).net_object
			ed = ucomplex(endderivativevalue).net_object
			z2 = _ComplexNumLib[_UncNumber].SplineIntegrate2(x2, y2, startboundary, sd, endboundary, ed)
			z = ucomplex(z2)
		else:
			sd = ufloat(startderivativevalue).net_object
			ed = ufloat(endderivativevalue).net_object
			z2 = _RealNumLib[_UncNumber].SplineIntegrate2(x2, y2, startboundary, sd, endboundary, ed)
			z = ufloat(z2)
		return z

	@staticmethod
	def fft(x):
		x2 = _asnetnarray(x, True)
		y2 = _ComplexNumLib[_UncNumber].Fft(x2)
		y = _fromnetnarray(y2)
		return y

	@staticmethod
	def ifft(x):
		x2 = _asnetnarray(x, True)
		y2 = _ComplexNumLib[_UncNumber].Ifft(x2)
		y = _fromnetnarray(y2)
		return y

	@staticmethod
	def dft(x):
		if _UncNumber is not _LinPropUncNumber:
			raise Exception("Dft supports only LinProp uncertainty objects")
		realTimeDomainData = not iscomplexarray(np.asarray(x))
		x2 = _asnetnarray(x, True)
		y2 = _ComplexUncNumLib[_UncNumber].Dft(x2, realTimeDomainData)
		y = _fromnetnarray(y2)
		return y

	@staticmethod
	def idft(x, realTimeDomainData=False):
		if _UncNumber is not _LinPropUncNumber:
			raise Exception("Idft supports only LinProp uncertainty objects")
		x2 = _asnetnarray(x, True)
		y2 = _ComplexUncNumLib[_UncNumber].Idft(x2, realTimeDomainData)
		y = _fromnetnarray(y2)
		if realTimeDomainData:
			y = umath.real(y)
		return y

	@staticmethod
	def roots(p):
		if _UncNumber is not _LinPropUncNumber:
			raise Exception("Roots supports only LinProp uncertainty objects")
		p2 = np.asarray(p)
		if p2.size < 2:
			d3 = np.array([])
		else:
			complexarray = iscomplexarray(p2)
			p3 = [_asnetnarray([[p2i]], True) for p2i in reversed(p2.flatten('F'))]
			res = _ComplexUncLinAlg[_UncNumber].NonLinearEig(p3)
			d2 = _fromnetnarray(res[1])
			if not complexarray and iszero(umath.imag(d2)):
				d3 = umath.real(d2)
			else:
				d3 = d2
		return d3

	@staticmethod
	def numerical_step(func, x, dx):
		if _UncNumber is not _LinPropUncNumber:
			raise Exception("Numerical step supports only LinProp uncertainty objects")
		x2 = _asnetnarray(x).data
		dx2 = _asnetdoublearray(dx)
		func2 = _NumericalFunctionDelegate(lambda _x: unumlib._numerical_function2(func, _x))
		y = _UncNumerical.NumericalStep(func2, x2, dx2, True)
		y2 = np.asarray([ufloat(yi) for yi in y])
		return y2

	@staticmethod
	def _numerical_function2(func, x):
		x2 = np.asarray([float(xi) for xi in x])
		y2 = func(x2)
		y = _asnetdoublearray(y2)
		return y

	@staticmethod
	def optimizer(func, xstart, p, covarianceweighting=True, weights=None, bndl=None, bndu=None, epsx=0.0, algorithm=OptimizerAlgorithm.LevenbergMarquardt):
		if _UncNumber is not _LinPropUncNumber:
			raise Exception("Optimizer supports only LinProp uncertainty objects")
		xstart2 = _asnetnarray(xstart).data
		p2 = _asnetnarray(p).data
		func2 = _ObjectiveFunction[_UncNumber](lambda _x, _p, _debug: unumlib._objective_function2(func, _x, _p, _debug))
		if bndl is None or bndu is None:
			if weights is None:
				x = _Optimizer.Start[_UncNumber](func2, xstart2, p2, epsx, covarianceweighting, algorithm)
			else:
				weights2 = _asnetdoublearray(weights)
				x = _Optimizer.StartWithWeights[_UncNumber](func2, xstart2, p2, weights2, epsx, covarianceweighting, algorithm)
		else:
			bndl2 = _asnetdoublearray(bndl)
			bndu2 = _asnetdoublearray(bndu)
			if weights is None:
				x = _Optimizer.StartWithBounds[_UncNumber](func2, xstart2, p2, bndl2, bndu2, epsx, covarianceweighting, algorithm)
			else:
				weights2 = _asnetdoublearray(weights)
				x = _Optimizer.StartWithWeightsAndBounds[_UncNumber](func2, xstart2, p2, weights2, bndl2, bndu2, epsx, covarianceweighting, algorithm)
		x2 = np.asarray([ufloat(xi) for xi in x])
		return x2

	@staticmethod
	def _objective_function2(func, x, p, debug):
		x2 = np.asarray([ufloat(xi) for xi in x])
		p2 = np.asarray([ufloat(pi) for pi in p])
		f2 = func(x2, p2)
		f = _asnetnarray(f2).data
		return f


class ustorage(object):
	@staticmethod
	def save_binary_file(x, filepath):
		x2 = _asnetobject(x)
		x2.BinarySerialize(filepath)

	@staticmethod
	def load_binary_file(filepath):
		try:
			x = _ComplexNArray[_UncNumber]().BinaryDeserialize(filepath)
		except:
			try:
				x = _RealNArray[_UncNumber]().BinaryDeserialize(filepath)
			except:
				try:
					x = _Complex[_UncNumber](0).BinaryDeserialize(filepath)
				except:
					try:
						x = _UncNumber(0).BinaryDeserialize(filepath)
					except:
						raise Exception("Wrong structure of binary file")
		x2 = _fromnetobject(x)
		return x2

	@staticmethod
	def save_xml_file(x, filepath):
		x2 = _asnetobject(x)
		x2.XmlSerialize(filepath)

	@staticmethod
	def load_xml_file(filepath):
		try:
			x = _ComplexNArray[_UncNumber]().XmlDeserialize(filepath)
		except:
			try:
				x = _RealNArray[_UncNumber]().XmlDeserialize(filepath)
			except:
				try:
					x = _Complex[_UncNumber](0).XmlDeserialize(filepath)
				except:
					try:
						x = _UncNumber(0).XmlDeserialize(filepath)
					except:
						raise Exception("Wrong structure of xml file")
		x2 = _fromnetobject(x)
		return x2

	@staticmethod
	def to_xml_string(x):
		x2 = _asnetobject(x)
		return x2.XmlSerializeToString()

	@staticmethod
	def from_xml_string(s):
		try:
			x = _ComplexNArray[_UncNumber]().XmlDeserializeFromString(s)
		except:
			try:
				x = _RealNArray[_UncNumber]().XmlDeserializeFromString(s)
			except:
				try:
					x = _Complex[_UncNumber](0).XmlDeserializeFromString(s)
				except:
					try:
						x = _UncNumber(0).XmlDeserializeFromString(s)
					except:
						raise Exception("Wrong structure of xml string")
		x2 = _fromnetobject(x)
		return x2

	@staticmethod
	def to_byte_array(x):
		x2 = _asnetobject(x)
		return bytes(x2.BinarySerializeToByteArray())

	@staticmethod
	def from_byte_array(b):
		b2 = _Array[_Byte](b)
		try:
			x = _ComplexNArray[_UncNumber]().BinaryDeserializeFromByteArray(b2)
		except:
			try:
				x = _RealNArray[_UncNumber]().BinaryDeserializeFromByteArray(b2)
			except:
				try:
					x = _Complex[_UncNumber](0).BinaryDeserializeFromByteArray(b2)
				except:
					try:
						x = _UncNumber(0).BinaryDeserializeFromByteArray(b2)
					except:
						raise Exception("Wrong structure of byte array")
		x2 = _fromnetobject(x)
		return x2


class uspecial(object):
	@staticmethod
	def linprop2mcprop(x):
		use_linprop()
		x_values = get_value(x)
		x_covariance = get_covariance(x)
		if isinstance(x, ufloat):
			use_mcprop()
			xmc = ufloat(x_values, np.sqrt(x_covariance[0, 0]))
		elif isinstance(x, ucomplex):
			use_mcprop()
			xmc = ucomplex(x_values, covariance=x_covariance)
		elif iscomplexarray(x):
			use_mcprop()
			xmc = ucomplexarray(x_values, x_covariance)
		else:
			use_mcprop()
			xmc = ufloatarray(x_values, x_covariance)
		return xmc
	
	@staticmethod
	def mcprop2linprop(ymc, xmc, x):
		use_mcprop()
		ymc2 = _fromnetnarray(_asnetnarray(_asunclist(ymc).data))
		xmc2 = _fromnetnarray(_asnetnarray(_asunclist(xmc).data))
		use_linprop()
		x2 = _fromnetnarray(_asnetnarray(_asunclist(x).data))
		xmc2sub = []
		x2sub = []
		for i in range(x2.size):
			if get_stdunc(x2[i]) > 0:
				xmc2sub.append(xmc2[i])
				x2sub.append(x2[i])
		use_mcprop()
		xmc_values = get_value(xmc2sub)
		xmc_covar = get_covariance(xmc2sub)
		ymc_values = get_value(ymc2)
		ymc_covar = get_covariance(ymc2)
		ymc_xmc_jacobi = get_jacobi2(ymc2, xmc2sub)
		use_linprop()
		x_covar = get_covariance(x2sub)
		xmc_x_jacobi = np.dot(np.linalg.cholesky(xmc_covar), np.linalg.inv(np.linalg.cholesky(x_covar)))
		x3 = [ufloatsystem(get_value(xmc_values[i]), x2sub, xmc_x_jacobi[i, :]) for i in range(len(x2sub))]
		y2 = [ufloatsystem(get_value(ymc_values[i]), x3, ymc_xmc_jacobi[i, :]) for i in range(ymc_values.size)]
		# Add non-linear contributions
		ymc_covar_lin = get_covariance(y2)
		ymc_covar_non_lin = ymc_covar - ymc_covar_lin
		if len(y2) == 1:
			l = ufloatarray(np.zeros(len(y2)), ymc_covar_non_lin, desc="Non-linear contribution")
		else:
			l = ufloatarray(np.zeros(len(y2)), ymc_covar_non_lin, desc="Non-linear contributions")
		y2 = y2 + l
		use_mcprop()
		if isinstance(ymc, ufloat):
			use_linprop()
			y = y2[0]
		elif isinstance(ymc, ucomplex):
			use_linprop()
			y = ucomplex(y2[0], imag=y2[1])
		elif iscomplexarray(ymc):
			use_linprop()
			y = [ucomplex(y2[2 * i], imag=y2[2 * i + 1]) for i in range(len(y2) // 2)]
			y = np.asarray(y).reshape(np.asarray(ymc).size, order='F')
		else:
			use_linprop()
			y = np.asarray(y2).reshape(np.asarray(ymc).size, order='F')
		return y


class ufloat(object):
	def __init__(self, value, stdunc=0.0, idof=0.0, id=None, desc=None):
		if isinstance(value, (int, float)):
			if stdunc != 0:
				id2, desc2 = _input_id_desc(id, desc)
				self._d = _UncNumber(float(value), float(stdunc), float(idof), id2, desc2)
			else:
				self._d = _UncNumber(float(value))
		elif type(value) is _UncNumber:
			self._d = value
		elif type(value) is _Real[_UncNumber]:
			self._d = value.Item
		elif isinstance(value, ufloat):
			self._d = value._d
		else:
			raise Exception("Unknown arguments")

	@staticmethod
	def _is_convertible(value):
		return (
			isinstance(value, (int, float)) |
			(type(value) is _UncNumber) |
			(type(value) is _Real[_UncNumber]) |
			isinstance(value, ufloat)
		)

	def __getstate__(self):
		state = ustorage.to_byte_array(self)
		return state

	def __setstate__(self, state):
		loaded = ustorage.from_byte_array(state)
		self._d = loaded._d

	def __repr__(self):
		return str(self.value) + _pm + str(self.stdunc)

	@property
	def value(self):
		return self._d.Value

	@property
	def fcn_value(self):
		return self._d.FcnValue

	@property
	def stdunc(self):
		return self._d.StdUnc

	@property
	def idof(self):
		return self._d.IDof

	@property
	def isconst(self):
		return self._d.IsConst

	@property
	def net_object(self):
		return self._d

	def get_moment(self, n):
		return self._d.GetMoment(int(n))

	def __pos__(self):
		return self

	def __neg__(self):
		return ufloat(self._d.Negative())

	def __add__(self, other):
		if iscomplex(other):
			return ucomplex(self) + other
		elif isinstance(other, np.ndarray):
			return np.asarray(self) + other
		elif ufloat._is_convertible(other):
			return ufloat(self._d.Add(ufloat(other)._d))
		else:
			return NotImplemented

	def __radd__(self, other):
		if iscomplex(other):
			return other + ucomplex(self)
		elif isinstance(other, np.ndarray):
			return other + np.asarray(self)
		elif ufloat._is_convertible(other):
			return ufloat(ufloat(other)._d.Add(self._d))
		else:
			return NotImplemented

	def __sub__(self, other):
		if iscomplex(other):
			return ucomplex(self) - other
		elif isinstance(other, np.ndarray):
			return np.asarray(self) - other
		elif ufloat._is_convertible(other):
			return ufloat(self._d.Subtract(ufloat(other)._d))
		else:
			return NotImplemented

	def __rsub__(self, other):
		if iscomplex(other):
			return other - ucomplex(self)
		elif isinstance(other, np.ndarray):
			return other - np.asarray(self)
		elif ufloat._is_convertible(other):
			return ufloat(ufloat(other)._d.Subtract(self._d))
		else:
			return NotImplemented

	def __mul__(self, other):
		if iscomplex(other):
			return ucomplex(self) * other
		elif isinstance(other, np.ndarray):
			return np.asarray(self) * other
		elif ufloat._is_convertible(other):
			return ufloat(self._d.Multiply(ufloat(other)._d))
		else:
			return NotImplemented

	def __rmul__(self, other):
		if iscomplex(other):
			return other * ucomplex(self)
		elif isinstance(other, np.ndarray):
			return other * np.asarray(self)
		elif ufloat._is_convertible(other):
			return ufloat(ufloat(other)._d.Multiply(self._d))
		else:
			return NotImplemented

	def __div__(self, other):
		if iscomplex(other):
			return ucomplex(self) / other
		elif isinstance(other, np.ndarray):
			return np.asarray(self) / other
		elif ufloat._is_convertible(other):
			return ufloat(self._d.Divide(ufloat(other)._d))
		else:
			return NotImplemented

	def __rdiv__(self, other):
		if iscomplex(other):
			return other / ucomplex(self)
		elif isinstance(other, np.ndarray):
			return other / np.asarray(self)
		elif ufloat._is_convertible(other):
			return ufloat(ufloat(other)._d.Divide(self._d))
		else:
			return NotImplemented

	def __truediv__(self, other):
		if iscomplex(other):
			return ucomplex(self) / other
		elif isinstance(other, np.ndarray):
			return np.asarray(self) / other
		elif ufloat._is_convertible(other):
			return ufloat(self._d.Divide(ufloat(other)._d))
		else:
			return NotImplemented

	def __rtruediv__(self, other):
		if iscomplex(other):
			return other / ucomplex(self)
		elif isinstance(other, np.ndarray):
			return other / np.asarray(self)
		elif ufloat._is_convertible(other):
			return ufloat(ufloat(other)._d.Divide(self._d))
		else:
			return NotImplemented

	def __pow__(self, other):
		if iscomplex(other):
			return ucomplex(self) ** other
		elif isinstance(other, np.ndarray):
			return np.asarray(self) ** other
		elif type(other) is int:
			return ufloat(self._d.Pow(other))
		elif get_fcn_value(self) < 0:
			return ucomplex(self) ** other
		elif ufloat._is_convertible(other):
			return ufloat(self._d.Pow(ufloat(other)._d))
		else:
			return NotImplemented

	def __rpow__(self, other):
		if iscomplex(other):
			return other ** ucomplex(self)
		elif isinstance(other, np.ndarray):
			return other ** np.asarray(self)
		elif ufloat._is_convertible(other):
			return ufloat(ufloat(other)._d.Pow(self._d))
		else:
			return NotImplemented

	def __abs__(self):
		return ufloat(self._d.Abs())

	def abs(self):
		return ufloat(self._d.Abs())

	def exp(self):
		return ufloat(self._d.Exp())

	def log(self):
		if get_fcn_value(self) < 0:
			return ucomplex(self).log()
		else:
			return ufloat(self._d.Log())

	def log10(self):
		if get_fcn_value(self) < 0:
			return ucomplex(self).log10()
		else:
			return ufloat(self._d.Log10())

	def pow(self, other):
		return self**other

	def sqrt(self):
		if get_fcn_value(self) < 0:
			return ucomplex(self).sqrt()
		else:
			return ufloat(self._d.Sqrt())

	def sin(self):
		return ufloat(self._d.Sin())

	def cos(self):
		return ufloat(self._d.Cos())

	def tan(self):
		return ufloat(self._d.Tan())

	def asin(self):
		return ufloat(self._d.Asin())

	def acos(self):
		return ufloat(self._d.Acos())

	def atan(self):
		return ufloat(self._d.Atan())

	def arcsin(self):
		return ufloat(self._d.Asin())

	def arccos(self):
		return ufloat(self._d.Acos())

	def arctan(self):
		return ufloat(self._d.Atan())

	def sinh(self):
		return ufloat(self._d.Sinh())

	def cosh(self):
		return ufloat(self._d.Cosh())

	def tanh(self):
		return ufloat(self._d.Tanh())

	def asinh(self):
		return ufloat(self._d.Asinh())

	def acosh(self):
		return ufloat(self._d.Acosh())

	def atanh(self):
		return ufloat(self._d.Atanh())

	def arcsinh(self):
		return ufloat(self._d.Asinh())

	def arccosh(self):
		return ufloat(self._d.Acosh())

	def arctanh(self):
		return ufloat(self._d.Atanh())

	def atan2(self, other):
		return ufloat(self._d.Atan2(ufloat(other)._d))

	def arctan2(self, other):
		return ufloat(self._d.Atan2(ufloat(other)._d))

	def sign(self):
		return self._d.Sign().Value

	def ellipk(self):
		return ufloat(self._d.Ellipk())

	def ellipe(self):
		return ufloat(self._d.Ellipe())

class ucomplex(object):
	def __init__(self, value, imag=0.0, covariance=None, idof=0.0, id=None, desc=None):
		if covariance is None:
			if iscomplex(value) and imag == 0:
				if isinstance(value, complex):
					self._d = _Complex[_UncNumber](_UncNumber(value.real), _UncNumber(value.imag))
				elif type(value) is _Complex[_UncNumber]:
					self._d = value
				else:
					self._d = value._d
			else:
				_real = ufloat(value)
				_imag = ufloat(imag)
				self._d = _Complex[_UncNumber](_real.net_object, _imag.net_object)
		else:
			_real = complex(value).real
			_imag = complex(value).imag + imag
			v = _Complex[_Number](_Number(_real), _Number(_imag))
			cv = _asnetnumbernarray(covariance)
			id2, desc2 = _input_id_desc(id, desc)
			self._d = _UncHelper.ComplexUncNumber(v, cv.Matrix, float(idof), id2, desc2)

	@staticmethod
	def _is_convertible(value):
		return (
			ufloat._is_convertible(value) | 
			iscomplex(value)
		)

	def __getstate__(self):
		state = ustorage.to_byte_array(self)
		return state

	def __setstate__(self, state):
		loaded = ustorage.from_byte_array(state)
		self._d = loaded._d

	def __repr__(self):
		return str(self.value) + _pm + str(self.stdunc)

	@property
	def real(self):
		return ufloat(self._d.Real())

	@property
	def imag(self):
		return ufloat(self._d.Imag())

	@property
	def value(self):
		return complex(self._d.Real().Value, self._d.Imag().Value)

	@property
	def fcn_value(self):
		return complex(self._d.Real().FcnValue, self._d.Imag().FcnValue)

	@property
	def stdunc(self):
		return complex(self._d.Real().StdUnc, self._d.Imag().StdUnc)
	
	@property
	def idof(self):
		return complex(self._d.Real().IDof, self._d.Imag().IDof)
	
	@property
	def isconst(self):
		return self._d.Real().IsConst and self._d.Imag().IsConst

	@property
	def net_object(self):
		return self._d

	def get_moment(self, n):
		return complex(self._d.Real().GetMoment(int(n)), self._d.Imag().GetMoment(int(n)))
	
	def __pos__(self):
		return self

	def __neg__(self):
		return ucomplex(self._d.Negative())

	def __add__(self, other):
		if isinstance(other, np.ndarray):
			return np.asarray(self) + other
		elif ucomplex._is_convertible(other):
			return ucomplex(self._d.Add(ucomplex(other)._d))
		else: 
			return NotImplemented

	def __radd__(self, other):
		if isinstance(other, np.ndarray):
			return other + np.asarray(self)
		elif ucomplex._is_convertible(other):
			return ucomplex(ucomplex(other)._d.Add(self._d))
		else: 
			return NotImplemented

	def __sub__(self, other):
		if isinstance(other, np.ndarray):
			return np.asarray(self) - other
		elif ucomplex._is_convertible(other):
			return ucomplex(self._d.Subtract(ucomplex(other)._d))
		else: 
			return NotImplemented

	def __rsub__(self, other):
		if isinstance(other, np.ndarray):
			return other - np.asarray(self)
		elif ucomplex._is_convertible(other):
			return ucomplex(ucomplex(other)._d.Subtract(self._d))
		else: 
			return NotImplemented

	def __mul__(self, other):
		if isinstance(other, np.ndarray):
			return np.asarray(self) * other
		elif ucomplex._is_convertible(other):
			return ucomplex(self._d.Multiply(ucomplex(other)._d))
		else: 
			return NotImplemented

	def __rmul__(self, other):
		if isinstance(other, np.ndarray):
			return other * np.asarray(self)
		elif ucomplex._is_convertible(other):
			return ucomplex(ucomplex(other)._d.Multiply(self._d))
		else: 
			return NotImplemented

	def __div__(self, other):
		if isinstance(other, np.ndarray):
			return np.asarray(self) / other
		elif ucomplex._is_convertible(other):
			return ucomplex(self._d.Divide(ucomplex(other)._d))
		else: 
			return NotImplemented

	def __rdiv__(self, other):
		if isinstance(other, np.ndarray):
			return other / np.asarray(self)
		elif ucomplex._is_convertible(other):
			return ucomplex(ucomplex(other)._d.Divide(self._d))
		else: 
			return NotImplemented

	def __truediv__(self, other):
		if isinstance(other, np.ndarray):
			return np.asarray(self) / other
		elif ucomplex._is_convertible(other):
			return ucomplex(self._d.Divide(ucomplex(other)._d))
		else: 
			return NotImplemented

	def __rtruediv__(self, other):
		if isinstance(other, np.ndarray):
			return other / np.asarray(self)
		elif ucomplex._is_convertible(other):
			return ucomplex(ucomplex(other)._d.Divide(self._d))
		else: 
			return NotImplemented

	def __pow__(self, other):
		if type(other) is int:
			return ucomplex(self._d.Pow(other))
		elif isinstance(other, np.ndarray):
			return np.asarray(self) ** other
		elif ucomplex._is_convertible(other):
			return ucomplex(self._d.Pow(ucomplex(other)._d))
		else: 
			return NotImplemented

	def __rpow__(self, other):
		if isinstance(other, np.ndarray):
			return other ** np.asarray(self)
		elif ucomplex._is_convertible(other):
			return ucomplex(ucomplex(other)._d.Pow(self._d))
		else: 
			return NotImplemented
	
	def __abs__(self):
		return ufloat(self._d.Abs())

	def abs(self):
		return ufloat(self._d.Abs())

	def exp(self):
		return ucomplex(self._d.Exp())

	def log(self):
		return ucomplex(self._d.Log())

	def log10(self):
		return ucomplex(self._d.Log10())

	def pow(self, other):
		return self**other

	def sqrt(self):
		return ucomplex(self._d.Sqrt())

	def sin(self):
		return ucomplex(self._d.Sin())

	def cos(self):
		return ucomplex(self._d.Cos())

	def tan(self):
		return ucomplex(self._d.Tan())

	def asin(self):
		return ucomplex(self._d.Asin())

	def acos(self):
		return ucomplex(self._d.Acos())

	def atan(self):
		return ucomplex(self._d.Atan())

	def arcsin(self):
		return ucomplex(self._d.Asin())

	def arccos(self):
		return ucomplex(self._d.Acos())

	def arctan(self):
		return ucomplex(self._d.Atan())

	def sinh(self):
		return ucomplex(self._d.Sinh())

	def cosh(self):
		return ucomplex(self._d.Cosh())

	def tanh(self):
		return ucomplex(self._d.Tanh())

	def asinh(self):
		return ucomplex(self._d.Asinh())

	def acosh(self):
		return ucomplex(self._d.Acosh())

	def atanh(self):
		return ucomplex(self._d.Atanh())

	def arcsinh(self):
		return ucomplex(self._d.Asinh())

	def arccosh(self):
		return ucomplex(self._d.Acosh())

	def arctanh(self):
		return ucomplex(self._d.Atanh())

	def conj(self):
		return ucomplex(self._d.Conj())

	def conjugate(self):
		return ucomplex(self._d.Conj())

	def angle(self, deg=False):
		if deg:
			return ufloat(self._d.Angle()) * (180.0 / np.pi)
		else:
			return ufloat(self._d.Angle())

	def phase(self, deg=False):
		return self.angle(deg)


class _const2014(object):
	"Physical Constants CODATA 2014"

	@property
	def _const(self):
		return _Const2014

	@property
	def deltavCs(self):
		"Hyperfine transition frequency of Cs-133 / Hz"
		return self._const.deltavCs

	@property
	def c0(self):
		"Speed of light in vacuum / (m/s)"
		return self._const.c0

	@property
	def mu0(self):
		"Vacuum magnetic permeability / (Vs/Am)"
		return self._const.mu0

	@property
	def ep0(self):
		"Vacuum electric permittivity / (As/Vm)"
		return self._const.ep0

	@property
	def Kcd(self):
		"Luminous efficacy / (lm/W)"
		return self._const.Kcd

	@property
	def Mu(self):
		"Molar mass _constant / (kg/mol)"
		return self._const.Mu

class _uconst2014(_const2014):
	"Physical Constants CODATA 2014"

	@property
	def _uconst(self):
		return _Const2014[_UncNumber]

	@property
	def G(self):
		"Newtonian constant of gravitation / (m^3/(kg*s^2))"
		return ufloat(self._uconst.G)

	@property
	def alpha(self):
		"Fine-structure constant"
		return ufloat(self._uconst.alpha)

	@property
	def Ryd(self):
		"Rydberg constant / (1/m)"
		return ufloat(self._uconst.Ryd)

	@property
	def mpsme(self):
		"Proton-electron mass ratio"
		return ufloat(self._uconst.mpsme)

	@property
	def Na(self):
		"Avogadro constant / (1/mol)"
		return ufloat(self._uconst.Na)

	@property
	def Kj(self):
		"Josephson constant / (Hz/V)"
		return ufloat(self._uconst.Kj)

	@property
	def k(self):
		"Boltzmann constant / (J/K)"
		return ufloat(self._uconst.k)

	@property
	def Rk(self): 
		"von Klitzing constant / Ohm"
		return ufloat(self._uconst.Rk)

	@property
	def e(self):
		"Elementary charge / C"
		return ufloat(self._uconst.e)

	@property
	def h(self):
		"Planck constant / Js"
		return ufloat(self._uconst.h)

	@property
	def me(self):
		"Electron mass / kg"
		return ufloat(self._uconst.me)

	@property
	def mp(self):
		"Proton mass / kg"
		return ufloat(self._uconst.mp)

	@property
	def u(self):
		"Atomic mass constant / kg"
		return ufloat(self._uconst.u)

	@property
	def F(self):
		"Faraday constant / (C/mol)"
		return ufloat(self._uconst.F)

	@property
	def R(self):
		"Molar gas constant / (J/(mol*K))"
		return ufloat(self._uconst.R)

	@property
	def eV(self):
		"Electron volt / J"
		return ufloat(self._uconst.eV)


class _const2014_90(_const2014):
	"Physical Constants CODATA 2014 for Conventional Electrical Units 90"

	@property
	def _const90(self):
		return _Const2014_90

	@property
	def Kj(self):
		"Conventional value of Josephson _constant / (Hz/V)"
		return self._const90.Kj

	@property
	def Rk(self):
		"Conventional value of von Klitzing _constant / Ohm"
		return self._const90.Rk

	@property
	def e(self):
		"Elementary charge / C"
		return self._const90.e

	@property
	def h(self):
		"Planck _constant / Js"
		return self._const90.h

class _uconst2014_90(_const2014_90):
	"Physical Constants CODATA 2014 for Conventional Electrical Units 90"

	@property
	def _uconst90(self):
		return _Const2014_90[_UncNumber]

	@property
	def Na(self):
		"Avogadro constant / (1/mol)"
		return ufloat(self._uconst90.Na)

	@property
	def F(self):
		"Faraday constant / (C/mol)"
		return ufloat(self._uconst90.F)

	@property
	def k(self):
		"Boltzmann constant / (J/K)"
		return ufloat(self._uconst90.k)


class _const2018(object):
	"Physical Constants CODATA 2018"

	@property
	def _const(self):
		return _Const2018

	@property
	def deltavCs(self):
		"Hyperfine transition frequency of Cs-133 / Hz"
		return self._const.deltavCs

	@property
	def c0(self):
		"Speed of light in vacuum / (m/s)"
		return self._const.c0

	@property
	def h(self):
		"Planck _constant / Js"
		return self._const.h

	@property
	def e(self):
		"Elementary charge / C"
		return self._const.e

	@property
	def k(self):
		"Boltzmann _constant / (J/K)"
		return self._const.k

	@property
	def Na(self):
		"Avogadro _constant / (1/mol)"
		return self._const.Na

	@property
	def Kcd(self):
		"Luminous efficacy / (lm/W)"
		return self._const.Kcd

	@property
	def Kj(self):
		"Josephson _constant / (Hz/V)"
		return self._const.Kj

	@property
	def Rk(self):
		"von Klitzing _constant / Ohm"
		return self._const.Rk

	@property
	def F(self):
		"Faraday _constant / (C/mol)"
		return self._const.F

	@property
	def R(self):
		"Molar gas _constant / (J/(mol*K))"
		return self._const.R

	@property
	def eV(self):
		"Electron volt / J"
		return self._const.eV

class _uconst2018(_const2018):
	"Physical Constants CODATA 2018"

	@property
	def _uconst(self):
		return _Const2018[_UncNumber]

	@property
	def G(self):
		"Newtonian _uconstant of gravitation / (m^3/(kg*s^2))"
		return ufloat(self._uconst.G)

	@property
	def alpha(self):
		"Fine-structure _uconstant"
		return ufloat(self._uconst.alpha)

	@property
	def mu0(self):
		"Vacuum magnetic permeability / (Vs/Am)"
		return ufloat(self._uconst.mu0)

	@property
	def ep0(self):
		"Vacuum electric permittivity / (As/Vm)"
		return ufloat(self._uconst.ep0)

	@property
	def Ryd(self):
		"Rydberg _uconstant / (1/m)"
		return ufloat(self._uconst.Ryd)

	@property
	def me(self):
		"Electron mass / kg"
		return ufloat(self._uconst.me)

	@property
	def are(self):
		"Electron relative atomic mass"
		return ufloat(self._uconst.are)

	@property
	def arp(self):
		"Proton relative atomic mass"
		return ufloat(self._uconst.arp)

	@property
	def mpsme(self):
		"Proton-electron mass ratio"
		return ufloat(self._uconst.mpsme)

	@property
	def mp(self):
		"Proton mass / kg"
		return ufloat(self._uconst.mp)

	@property
	def u(self):
		"Atomic mass _uconstant / kg"
		return ufloat(self._uconst.u)

	@property
	def Mu(self):
		"Molar mass _uconstant / (kg/mol)"
		return ufloat(self._uconst.Mu)


uconst2014 = _uconst2014()
"Physical Constants CODATA 2014"

uconst2014_90 = _uconst2014_90()
"Physical Constants CODATA 2014 for Conventional Electrical Units 90"

uconst2018 = _uconst2018()
"Physical Constants CODATA 2018"

uconst = _uconst2018()
"Newest Physical Constants"
