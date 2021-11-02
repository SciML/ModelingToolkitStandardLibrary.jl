# ModelingToolkitStandardLibrary.jl

A standard library of components to model the world and beyond.

- [Electrical modeling](#electrical-modeling)
  - [Basic electrical components](#basic-electrical-components)
    - [Ground](#ground)
    - [Resistor](#resistor)
    - [Capacitor](#capacitor)
  - [Ideal electrical elements](#ideal-electrical-elements)
    - [Short](#short)
    - [IdealOpAmp](#idealopamp)
  - [Sensors](#sensors)
    - [CurrentSensor](#currentsensor)
    - [PotentialSensor](#potentialsensor)
    - [VoltageSensor](#voltagesensor)
  - [Voltage/current sources](#voltagecurrent-sources)
    - [ConstantVoltage](#constantvoltage)
    - [StepVoltage](#stepvoltage)
    - [SineVoltage](#sinevoltage)
- [Thermal modeling](#thermal-modeling)
  - [Basic thermal components](#basic-thermal-components)
    - [Thermal ground](#thermal-ground)
    - [Heat capacitor](#heat-capacitor)
    - [Thermal conductor](#thermal-conductor)
    - [Thermal resistor](#thermal-resistor)
    - [Convective conductor](#convective-conductor)
    - [Convective resistor](#convective-resistor)
    - [Thermal radiation](#thermal-radiation)
    - [Thermal collector](#thermal-collector)
  - [Thermal sensors](#thermal-sensors)
    - [Temperature sensor](#temperature-sensor)
    - [Relative temperature sensor](#relative-temperature-sensor)
    - [Heat flow sensor](#heat-flow-sensor)
  - [Thermal sources](#thermal-sources)
    - [Fixed heat flow](#fixed-heat-flow)
    - [Fixed temperature](#fixed-temperature)

## Electrical modeling

Currently, ModelingToolkitStandardLibrary contains basic electrical components, ideal electrical elements, sensors, and voltage/current sources.

### Basic electrical components

#### Ground

**Function**: `Ground(;name)`

**Description**: Ground node with the potential of zero and connector `g`.
Note that specifying the macro `@named sys = Ground()` is equivalent to setting
`sys = Ground(;name, sys)`. Either method will suffice, and there is no need to
type the name twice. The same principle applies to the other electrical components.

**Connectors**:
- `g`

### Resistor

**Function**: `Resistor(;name, R = 1.0)`

**Observables**:
- `R`: resistance (negative, zero, positive)

**States**
- `v`: voltage

**Connectors**:
- positive pin
- negative pin

### Capacitor

**Function**: `Capacitor(;name, C = 1.0)`

**Observables**:
- `C`: capacitance (zero or positive)

**Connectors**:
- positive pin
- negative pin

---

### Ideal electrical elements

#### Short

**Function**: `Short()`

**Description**: Short cut branch.

**Connectors**:
- positive pin
- negative pin

#### IdealOpAmp

**Function**: `IdealOpAmp(;name)`

**Description**: The ideal operational amplifier.

**States**:
- `v1`: voltage of the left port
- `v2`: voltage of the right port
- `i1`: current of the left port
- `i2`: current of the right port

**Connectors**:
- positive pin (left port)
- negative pin (left port)
- positive pin (right port)
- negative pin (right port)

---

### Sensors

#### CurrentSensor

**Function**: `CurrentSensor(;name)`

**States**
- `i`: current value from the positive to the negative pin

**Connectors**:
- positive pin
- negative pin


#### PotentialSensor

**Function**: `PotentialSensor(;name)`

**Connectors**:
- pin (which is to be measured)

#### VoltageSensor

**Function**: `VoltageSensor(;name)`

**States**:
- `v`: value of voltage between the two pins

**Connectors**:
- positive pin
- negative pin

---

### Voltage/current sources

#### ConstantVoltage

**Function**: `ConstantVoltage(;name, V = 1.0)`

**Description**: The source for an ideal constant voltage.

**Observables**:
- `V`: value of constant voltage

**Connectors**:
- positive pin
- negative pin

#### StepVoltage

**Function**: `StepVoltage(;name, offset = 0.0, starttime = 0.0, height = 0.0)`

**Description**: Step voltage source.

**Observables**:
- `offset`: voltage offset
- `starttime`: time offset
- `height`: height of the step

**Connectors**:
- positive pin
- negative pin

#### SineVoltage

**Function**: `SineVoltage(;name, offset = 0.0, amplitude = 0.0, frequency = 0.0, starttime = 0.0, phase = 0.0)`

**Description**: Sine voltage source.

**Observables**:
- `offset`: voltage offset
- `amplitude`: amplitude of the sine wave
- `frequency`: frequency of the sine wave
- `starttime`: time offset
- `phase`: phase of the sine wave

**Connectors**:
- positive pin
- negative pin

## Thermal modeling

ModelingToolkitStandardLibrary contains basic thermal components for modeling heat transfer and fluid heat flow.

### Basic thermal components

#### Thermal ground

**Function**: `ThermalGround(;name)`

**Description**: Thermal port for 1-dimensional heat transfer with temperature set to zero.
Note that specifying the macro `@named sys = ThermalGround()` is equivalent to setting
`sys = ThermalGround(;name, sys)`. Either method will suffice, and there is no need to
type the name twice. The same principle applies to the other thermal components.

#### Heat capacitor

**Function**: `HeatCapacitor(;name, C = 1.0)`

**Observables**:
- `C`: heat capacity (zero or positive)

**State**:
- `T`: temperature (in Kelvin)

**Connectors**:
- heat port

#### Thermal conductor

**Function**: `ThermalConductor(;name, G = 1.0)`

**Observables**:
- `G`: thermal conductance

**States**:
- `Q_flow`: heat flow rate
- `T`: temperature (in Kelvin)

**Connectors**:
- two heat ports

#### Thermal resistor

**Function**: `ThermalResistor(;name, R = 1.0)`

**Description**: The model operates on the same principle as `ThermalConductor`, but relies on thermal resistance instead of thermal conductance.

**Observables**:
- `R`: thermal resistance

**States**:
- `Q_flow`: heat flow rate
- `T`: temperature (in Kelvin)

**Connectors**:
- two heat ports

#### Convective conductor

**Function**: `ConvectiveConductor(;name; G = 1.0)`

**Description**: Model of linear heat convection.

**Observables**:
- `G`: convective thermal conductance

**States**:
- `Q_flow`: heat flow rate
- `dT`: temperature difference (in Kelvin)

**Connectors**:
- two heat ports (for modeling of the fluid flow over the solid)

#### Convective resistor

**Function**: `ConvectiveResistor(;name; R = 1.0)`

**Description**: Model of linear heat convection. Works analogously to the above model, but relies on convective thermal resistance instead of convective thermal conductance.

**Observables**:
- `R`: convective thermal resistance

**States**:
- `Q_flow`: heat flow rate
- `dT`: temperature difference (in Kelvin)

**Connectors**:
- two heat ports (for modeling of the fluid flow over the solid)

#### Thermal radiation

**Function**: `BodyRadiation(;name; G = 1.0)`

**Description**: Thermal radiation model.

**Observables**:
- `G`: net radiation conductance between two surfaces
- Stefan-Boltzmann constant

**States**:
- `Q_flow`: heat flow rate

**Connectors**:
- two heat ports

#### Thermal collector

**Function**: `ThermalCollector(;name, N = 1)`

**Description**: A model for collecting the heat flows from multiple heatports
to a singular heatport.

**Observables**:
- `N`: the number of input heatports

**States**:
- `Q_flow`: heat flow rate
- `T`: temperature (in Kelvin)

**Connectors**:
- `hp1...hpN`: the respective heatports
- `collector_port`: the target collector heatport

---

### Thermal sensors

#### Temperature sensor

**Function**: `TemperatureSensor(;name)`

**Description**: Ideal absolute temperature sensor which outputs the temperature (in Kelvin) of the connected port.

**States**:
- `T`: temperature (in Kelvin)

**Connectors**:
- heat port

#### Relative temperature sensor

**Function**: `RelativeTemperatureSensor(;name)`

**Description**: The output of the sensor is the relative temperature, i.e., the difference of the two ports, given in Kelvin.

**States**:
- `T`: temperature (in Kelvin)

**Connectors**:
- two heat ports

#### Heat flow sensor

**Function**: `HeatFlowSensor(;name)`

**Description**: The model monitors the heat flow rate of the component. Its output is positive when the direction the heat flow is from the first port to the second one.

**States**:
- `Q_flow`: heat flow rate

**Connectors**:
- two heat ports

---

### Thermal sources

#### Fixed heat flow

**Function**: `FixedHeatFlow(;name, Q_flow=1.0, T₀=293.15, α=0.0)`

**Observables**:
- `Q_flow`: the constant amount of heat flow rate
- `T₀`: the reference temperature
- `α`: this parameter simulates temperature-dependent loss (if the specified value is not 0)

**Connectors**:
- heat port

#### Fixed temperature

**Function**: `FixedTemperature(;name, T = 0.0)`

**Description**: The model defines a fixed temperature (in Kelvin) at a given port.

**Observables**:
- `T`: temperature (in Kelvin)

**Connectors**:
- heat port
