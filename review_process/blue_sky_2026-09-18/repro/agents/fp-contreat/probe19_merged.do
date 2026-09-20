* driver: runs probes 16, 17, 18 in one Stata instance (absolute paths on purpose:
* each sub-file issues `clear all')
cap noisily do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat/probe16_poly_lat_uz.do"
dis "SECTION16_RC=" _rc
cap noisily do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat/probe17_within_collapse.do"
dis "SECTION17_RC=" _rc
cap noisily do "/private/tmp/claude-501/-Users-hegland-Documents-Research-wooldid-dev/f98806d6-cb0a-446e-9291-a575a3bcab42/scratchpad/bluesky/probes/fp-contreat/probe18_unbalanced.do"
dis "SECTION18_RC=" _rc
dis "ALLDONE"
