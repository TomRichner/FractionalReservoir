# fig_ied_stim_blocks statistics

Run: ied_stim_blocks, amplitude 0.2, 12 networks (seeds [42 73 101 137 211 307 409 503 601 709 811 907]), 4 no-stimulation + 4 stimulation blocks per network, 300 s analysed per block after 10 s settle (plus the model warmup); noise seed per block base.noise_seed + (seed - 42) + 1000*block.

Example traces: no stimulation seed 42 block 1, 225-255 s of the analysed block; stimulation seed 42 block 5, 50-80 s of the analysed block.

## Candidate events per minute, pooled per condition, and normalised block rates

| Network | Seed | No stim (ev/min) | Stim (ev/min) | Stim / no stim | No-stim blocks (% of pooled) | Stim blocks (% of pooled no stim) |
|---|---|---|---|---|---|---|
| N1 | 42 | 18.4 | 2.1 | 0.11 | 108, 107, 95, 89 | 10, 13, 10, 13 |
| N2 | 73 | 9.8 | 1.4 | 0.15 | 92, 109, 94, 105 | 16, 14, 12, 16 |
| N3 | 101 | 23.5 | 1.9 | 0.08 | 93, 107, 101, 99 | 6, 7, 10, 9 |
| N4 | 137 | 13.3 | 1.4 | 0.10 | 99, 105, 99, 96 | 8, 5, 18, 11 |
| N5 | 211 | 12.8 | 2.4 | 0.18 | 86, 106, 106, 102 | 17, 19, 17, 20 |
| N6 | 307 | 39.6 | 9.1 | 0.23 | 101, 101, 102, 95 | 21, 25, 24, 22 |
| N7 | 409 | 21.9 | 2.5 | 0.12 | 98, 106, 102, 94 | 13, 13, 10, 11 |
| N8 | 503 | 25.4 | 5.1 | 0.20 | 97, 99, 108, 95 | 19, 16, 23, 23 |
| N9 | 601 | 21.4 | 3.1 | 0.15 | 104, 96, 111, 88 | 20, 10, 11, 18 |
| N10 | 709 | 54.1 | 10.3 | 0.19 | 101, 101, 98, 101 | 20, 18, 21, 18 |
| N11 | 811 | 28.9 | 10.4 | 0.36 | 101, 101, 89, 110 | 37, 33, 35, 40 |
| N12 | 907 | 15.3 | 1.0 | 0.07 | 102, 102, 102, 95 | 4, 4, 12, 7 |

Ratio stimulation / no stimulation: median 0.15, range [0.07, 0.36]; lower under stimulation in 12 of 12 networks.
Events per block: no stimulation median 108 [45, 273]; stimulation median 12 [3, 57]; zero-event blocks: 0 no stim, 0 stim.

## State change with stimulation (%; per network from block means), median [min, max]

- Mean rate: +53.4 [+41.7, +64.1]
- SFA feedback: +53.4 [+41.7, +64.0]
- Depression 1 - b_1: +24.1 [+16.0, +31.8]
