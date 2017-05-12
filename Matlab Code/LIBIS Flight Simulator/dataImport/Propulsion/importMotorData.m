function motorData = importMotorData()

% AXI 5330/18 GOLD LINE V2
    motorData(1).Model = 'AXI 5330/18 GOLD LINE V2';
    motorData(1).I0    = 1.9;   %A (No load current)
    motorData(1).Kv    = 259;   %RPM/V
    motorData(1).Rm    = 32e-3; %Ohm (Internal resistance)
    motorData(1).MaxEfficiencyCurrent = [25, 60]; %Intensity range for eta>85%
    motorData(1).CurrentCapacity = 76; %Max. intensity (A) in 60s
    motorData(1).MaxEfficiency = 0.92;
    motorData(1).MaxPower = 2870; %W
    motorData(1).Weight = 0.672; %kg


% AXI 5320/28 GOLD LINE V2
    motorData(2).Model = 'AXI 5320/28 GOLD LINE V2';
    motorData(2).I0   = 1.3;   %A (No load current)
    motorData(2).Kv   = 249;   %RPM/V
    motorData(2).Rm   = 57e-3; %Ohm (Internal resistance)
    motorData(2).MaxEfficiencyCurrent = [10, 36]; %Intensity range for eta>85%
    motorData(2).CurrentCapacity = 76; %Max. intensity (A) in 60s
    motorData(2).MaxEfficiency = 0.94;
    motorData(2).MaxPower = 1600; %W
    motorData(2).Weight = 0.515; %kg
    
end