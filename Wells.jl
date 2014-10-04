module Wells

using Linv

stehfestcoefficients = Linv.getstehfestcoefficients()

function makedrawdownwithconstantheadboundary(drawdown::Function) # constant head boundary
	return (R::Number, args...)->begin
		t = args[1]
		r = args[2]
		otherargs = args[3:end]
		return drawdown(t, r, otherargs...) - drawdown(t, R - r, otherargs...)
	end
end

function makedrawdownwithzerofluxboundary(drawdown::Function) # zero flux boundary
	return (R::Number, args...)->begin
		t = args[1]
		r = args[2]
		otherargs = args[3:end]
		return drawdown(t, r, otherargs...) + drawdown(t, R - r, otherargs...)
	end
end

function Ei(u::Number) # Thies function
	if( u > 60 || u <= 0 )
		return 0.0
	end
	if u < 1
		s1 = ((((0.00107857 * u - 0.00976004) * u + 0.05519968) * u - 0.24991055) * u + 0.99999193) * u - 0.57721566;
		return s1 - log(u);
	else
		s1 = (((u + 8.57332875) * u + 18.05901697) * u + 8.63476088) * u + 0.267773734;
		s2 = (((u + 9.57332875) * u + 25.63295615) * u + 21.09965308) * u + 3.958496923;
		return (1 / u) * (s1 / s2) * exp(-u);
	end
end

# Hantush Function W
# calculated by logarithmic integration using Simpson's rule
# source: W.Kinzelbach - Groundwater modelling
function Whantush( ra::Number, rb::Number )
        if( ra < 1e-20 )
                return 2 * K0( rb )
        end
        if( ra > 2000.0 )
                return 0
        end
        ug = log( ra )
        og = 10.0
        hi = ( og - ug ) / 24
        sub1 = rb * rb / 4.0
        sub2 = hi * 2
        x4 = ug + hi
        x2 = ug
        s4 = s2 = 0
        i = 0
        while true
                s4 += exp( - exp( x4 ) - sub1 / exp( x4 ) )
                x4 += sub2
                if( i == 10 )
                        break
                end
                x2 += sub2
                s2 += exp( - exp( x2 ) - sub1 / exp( x2 ) )
                i++
        end
        return hi * ( exp( - exp( ug ) - sub1 / exp( ug ) ) + 4 * s4 + 2 * s2 + exp( - exp( og ) - sub1 / exp( og ) ) ) / 3
end

function theisdrawdown(t::Number, r::Number, T::Number, S::Number, Q::Number) # constant pumping rate
	u = r ^ 2 * S / (4 * T * t)
	return Q * Ei(u) / (4 * pi * T)
end

function theisdrawdown(t::Number, r::Number, T::Number, S::Number, Qm::Matrix) # step-wise changes in the pumping rate
	dd = 0.
	Qprev = 0.
	Qtime = Qm[1:end, 1] # first time
	Q = Qm[1:end, 2] # next pumping rate
	# println( "time ", Qtime )
	# println( "rate ", Q )
	i = 1
	while i <= size(Qtime)[1] && t > Qtime[i]
		dd += theisdrawdown(t - Qtime[i], r, T, S, Q[i] - Qprev)
		Qprev = Q[i]
		i += 1
	end
	return dd
end

#TODO we need to add the laplace solution for any functional form of the pumping rate (currently in well.c)

function runtheistests()
	T = 100
	S = 0.02
	Q = 2
	r = 10
	t = 3600:3600*24*365:3600*24*365*10
	for i in 1:size(t)[1]
		println(theisdrawdown(t[i], r, T, S, Q))
	end
	theisdrawdownwithzerofluxboundary = makedrawdownwithzerofluxboundary(theisdrawdown)
	R = map(i->i * r, 2:5)
	for i in 1:size(t)[1]
		println(map(R->theisdrawdownwithzerofluxboundary(R, t[i], r, T, S, Q), R))
	end
end

# runtheistests()

function K0(x::Number) # zerorh order bessel function of second kind
	return besselk(0, x)
end

# avci parameters -- note that 1 refers to the upper aquifer and 2 refers to the lower aquifer
# K -- Permeability
# L -- aquifer thickness
# Sc -- aquifer specific storage coefficient
# ra -- radius of the leaky well
# R -- distance from injection well to leaky well
# omega -- leaky well resistivity
# deltah = h1 - h2 where h1 is the head in the upper aquifer and h2 is the head in the lower aquifer
# r1 -- distance from leaky well to observation well in the upper aquifer
# r2 -- distance from leaky well to observation well in the lower aquifer
function laplaceavciflow(Qw::Number, K1::Number, K2::Number, L1::Number, L2::Number, Sc1::Number, Sc2::Number, ra::Number, R::Number, omega::Number, deltah::Number, s::Number)
	numerator = Qw * K0(R * sqrt(s * Sc2 / K2)) / (2 * pi * K2 * L2) + deltah
	denominator = s * (omega + K0(ra * sqrt(s * Sc2 / K2)) / (2 * pi * K2 * L2) + K0(ra * sqrt(s * Sc1 / K1)) / (2 * pi * K1 * L1))
	return numerator / denominator
end

function makeavciflow(Qw::Number, K1::Number, K2::Number, L1::Number, L2::Number, Sc1::Number, Sc2::Number, ra::Number, R::Number, omega::Number, deltah::Number)
	return Linv.makelaplaceinverse(s->laplaceavciflow(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, s))
end

function avciflow(Qw::Number, K1::Number, K2::Number, L1::Number, L2::Number, Sc1::Number, Sc2::Number, ra::Number, R::Number, omega::Number, deltah::Number, t::Number)
	return Linv.linv(s->laplaceavciflow(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, s), stehfestcoefficients, t)
end

# the next three functions give the head buildup in the upper aquifer
function laplaceavcideltahead1(Qw::Number, K1::Number, K2::Number, L1::Number, L2::Number, Sc1::Number, Sc2::Number, ra::Number, R::Number, omega::Number, deltah::Number, r1::Number, s::Number)
	return laplaceavciflow(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, s) * K0(r1 * sqrt(s * Sc1 / K1)) / (2 * pi * K1 * L1)
end

function makeavcideltahead1(Qw::Number, K1::Number, K2::Number, L1::Number, L2::Number, Sc1::Number, Sc2::Number, ra::Number, R::Number, omega::Number, deltah::Number, r1::Number)
	return Linv.makelaplaceinverse(s->laplaceavcideltahead1(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, r1, s))
end

function avcideltahead1(Qw::Number, K1::Number, K2::Number, L1::Number, L2::Number, Sc1::Number, Sc2::Number, ra::Number, R::Number, omega::Number, deltah::Number, r1::Number, t::Number)
	return Linv.linv(s->laplaceavcideltahead1(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, r1, s), stehfestcoefficients, t)
end

# the next three functions give the head buildup in the lower aquifer
function laplaceavcideltahead2(Qw::Number, K1::Number, K2::Number, L1::Number, L2::Number, Sc1::Number, Sc2::Number, ra::Number, R::Number, omega::Number, deltah::Number, r2::Number, rw::Number, s::Number)
	return -laplaceavciflow(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, s) * K0(r2 * sqrt(s * Sc2 / K2)) / (2 * pi * K2 * L2) + Qw * K0(rw * sqrt(s * Sc2 / K2)) / (2 * pi * K2 * L2 * s)
end

function makeavcideltahead2(Qw::Number, K1::Number, K2::Number, L1::Number, L2::Number, Sc1::Number, Sc2::Number, ra::Number, R::Number, omega::Number, deltah::Number, r2::Number, rw::Number)
	return Linv.makelaplaceinverse(s->laplaceavcideltahead2(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, r2, rw, s))
end

function avcideltahead2(Qw::Number, K1::Number, K2::Number, L1::Number, L2::Number, Sc1::Number, Sc2::Number, ra::Number, R::Number, omega::Number, deltah::Number, r2::Number, rw::Number, t::Number)
	return Linv.linv(s->laplaceavcideltahead2(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, r2, rw, s), stehfestcoefficients, t)
end

function runavcitests()
	Qw = .1 # m^3/sec
	K1 = 1e-3 # m/sec -- pervious
	K2 = 1e-5 # m/sec -- semi-pervious
	L1 = 100 # m
	L2 = 200 # m
	Sc1 = 7e-5 # m^-1 -- dense, sandy gravel
	Sc2 = 1e-5 # m^-1 -- fissured rock
	ra = .1 # m
	R = 100 # m
	omega = 1e3 # no resistance
	deltah = 0 # m
	r1 = 50 # m
	r2 = 100 # m
	rw = 25 # m

	af = makeavciflow(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah)
	adh1 = makeavcideltahead1(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, r1)
	adh2 = makeavcideltahead2(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, r2, rw)

	# println(map(af, 3600:3600:3600*24))
	t = 3600:3600*24*365:3600*24*365*10 # time in seconds; range from 1 hour to 10 years for every year
	vals = map(af, t)
	dh1s = map(adh1, t)
	dh2s = map(adh2, t)
	for i in 1:size(t)[1]
		time = t[i]
		v = vals[i]
		dh1 = dh1s[i]
		dh2 = dh2s[i]
		println("$time: $v, $dh1, $dh2")
	end
end

# runavcitests()

end
