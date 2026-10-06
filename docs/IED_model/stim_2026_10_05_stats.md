# fig_ied_stim statistics

Selected amplitude 0.2; 8 realizations (seeds [42 73 101 137 211 307 409 503]).

## Candidate events per minute by input level (median [min, max])

| Level | Per-window detector | Fixed criterion |
|---|---|---|
| 0 | 18.6 [9.6, 36.0] | 18.6 [9.6, 36.0] |
| 0.05 | 13.2 [2.4, 26.4] | 10.8 [1.2, 21.6] |
| 0.1 | 7.2 [1.2, 18.0] | 5.4 [0.0, 18.0] |
| 0.2 | 2.4 [0.0, 10.8] | 0.6 [0.0, 16.8] |
| 0.3 | 0.0 [0.0, 3.6] | 0.0 [0.0, 9.6] |

## At the selected amplitude (median [min, max] across realizations; event-shape rows omit realizations with no events)

| Quantity | No stimulation | Stimulation | Change |
|---|---|---|---|
| Events/min | 18.6 [9.6, 36] | 2.4 [0, 10.8] | -15 |
| Mean rate | 0.129 [0.0935, 0.163] | 0.195 [0.154, 0.232] | +0.0676 |
| SFA feedback | 0.0644 [0.0468, 0.0817] | 0.0977 [0.0769, 0.116] | +0.0338 |
| Synaptic resources | 0.614 [0.535, 0.707] | 0.524 [0.459, 0.613] | -0.0956 |
| Within-subnet correlation | 0.249 [0.161, 0.374] | 0.158 [0.125, 0.185] | -0.0926 |
| Between-subnet correlation | 0.00446 [-0.00606, 0.027] | -0.000122 [-0.0044, 0.00527] | -0.00415 |
| Median recruitment | 0.0921 [0.0814, 0.14] | 0.0714 [0.06, 0.131] | -0.0171 |
| Waveform FWHM (s) | 0.239 [0.192, 0.325] | 0.22 [0.178, 0.265] | -0.025 |

Event-rate ratio stimulation / no stimulation: median 0.12, range [0.00, 0.30]; below 1 in 8 of 8.

## Band power change (dB), median [min, max]

- Delta: -3.91 [-6.97, -0.89]; below 0 in 8 of 8
- Theta: +0.07 [-1.86, +0.87]; below 0 in 4 of 8
- Alpha: +0.35 [-0.14, +0.63]; below 0 in 1 of 8
- Beta: +0.12 [-0.11, +0.22]; below 0 in 3 of 8

## State change (%), median [min, max]

- Mean rate: +54.3 [+42.1, +64.4]
- SFA feedback: +54.3 [+42.1, +64.3]
- Synaptic resources: -14.6 [-16.6, -12.7]
