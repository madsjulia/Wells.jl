import Wells
using Base.Test

function runmonotonicitytest()
	T = 100
	S = 0.02
	Q = 2
	r = 10
	R = 100
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
	lambda = 100 * r
	ts = 3600:3600*24:3600*24*365*10
	fs = Function[]
	push!(fs, t->Wells.theisdrawdown(t, r, T, S, Q))
	theisdrawdownwithzerofluxboundary = Wells.makedrawdownwithzerofluxboundary(Wells.theisdrawdown)
	push!(fs, t->theisdrawdownwithzerofluxboundary(R, t, r, T, S, Q))
	push!(fs, Wells.makeavcideltahead1(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, r1))
	push!(fs, Wells.makeavcideltahead2(Qw, K1, K2, L1, L2, Sc1, Sc2, ra, R, omega, deltah, r2, rw))
	push!(fs, t->Wells.hantushleakydrawdown(t, r, T, S, Q, lambda))
	lastdrawdowns = zeros(length(fs))
	for t in ts
		for j = 1:length(fs)
			thisdrawdown = fs[j](t)
			@test thisdrawdown >= lastdrawdowns[j]
			lastdrawdowns[j] = thisdrawdown
		end
	end
end

function hantushlimittest()
	T = 100
	S = 0.02
	Q = 2
	r = 10
	lambda = 1e6#theis and leakyhantush should be the same for large lambda
	ts = 0:3600*24:3600*24*365*10
	for t in ts
		@test_approx_eq_eps Wells.theisdrawdown(t, r, T, S, Q) Wells.hantushleakydrawdown(t, r, T, S, Q, lambda) 1e-6
	end
end

function timedepmacrotest()
	T = 100
	S = 0.02
	Q = 2
	r = 10
	Qm = Array(Float64, T, 2)
	for i = 1:size(Qm, 1)
		Qm[i, 1] = i#set the time
		Qm[i, 2] = 1 + 0.5 * randn()
	end
	for t in linspace(0, T, 2 * T)
		@test Wells.theisdrawdown(t, r, T, S, Qm) == Wells.theisdrawdownmanual(t, r, T, S, Qm)
	end
	@time for t in linspace(0, T, 2 * T * 1000)
		Wells.theisdrawdown(t, r, T, S, Qm)
	end
	@time for t in linspace(0, T, 2 * T * 1000)
		Wells.theisdrawdownmanual(t, r, T, S, Qm)
	end
end

runmonotonicitytest()
hantushlimittest()
timedepmacrotest()
