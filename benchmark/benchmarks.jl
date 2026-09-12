using BenchmarkTools
using ModelingToolkit
using ModelingToolkit: t_nounits as t
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Blocks: Sine
using OrdinaryDiffEq

const SUITE = BenchmarkGroup()

function rc_circuit()
    @named source = Sine(offset = 1, amplitude = 10, frequency = 5)
    @named voltage = Voltage()
    @named resistor = Resistor(R = 1)
    @named capacitor = Capacitor(C = 1, v = 0.0)
    @named ground = Ground()
    @named voltage_sensor = VoltageSensor()
    @named current_sensor = CurrentSensor()

    connections = [
        connect(source.output, voltage.V)
        connect(voltage.p, resistor.p)
        connect(resistor.n, current_sensor.p)
        connect(current_sensor.n, capacitor.p)
        connect(capacitor.n, voltage.n, ground.g)
        connect(capacitor.p, voltage_sensor.p)
        connect(capacitor.n, voltage_sensor.n)
    ]

    return System(
        connections, t;
        systems = [
            resistor, capacitor, source, voltage, ground,
            voltage_sensor, current_sensor,
        ], name = :rc_circuit
    )
end

SUITE["rc_circuit"] = BenchmarkGroup()
SUITE["rc_circuit"]["build"] = @benchmarkable rc_circuit()
model = rc_circuit()
SUITE["rc_circuit"]["mtkcompile"] = @benchmarkable mtkcompile($model)
sys = mtkcompile(model)
SUITE["rc_circuit"]["odeproblem"] = @benchmarkable ODEProblem($sys, [], (0.0, 10.0))
prob = ODEProblem(sys, [], (0.0, 10.0))
SUITE["rc_circuit"]["solve"] = @benchmarkable solve(
    $prob, Tsit5(); save_everystep = false, abstol = 1.0e-8, reltol = 1.0e-8
)
