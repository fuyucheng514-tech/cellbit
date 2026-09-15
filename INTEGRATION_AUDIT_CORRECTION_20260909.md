# Audit correction note

Date (UTC): 2026-09-09T14:07:02Z
Original audit: /home/data/fyc/dna2bit_sag_original_integrated_20260909/INTEGRATION_AUDIT_20260909_ATTEMPT_002.md
Original audit SHA256: f86662f147eb2e9a1226c0b638a3bf9c2eb4bc2aec7b4240d9c0a960aac275a9

The sentence in the original attempt-002 audit saying that attempts 001, 002, 003 and 004 were failed is corrected here. The accurate status is:

- smoke attempt_001: failed receipt/path assertion, preserved;
- smoke attempt_002: failed top-level sketch-receipt assertion, preserved;
- smoke attempt_003: successful Stage1/3A main smoke and packed-search comparison;
- smoke attempt_004: failed Stage3B smoke because its temporary test script still used the stale CheckM2 path, preserved;
- smoke attempt_005: successful Stage3B upstream, Stage3B tractor, final view and final CheckM2 smoke.

No scientific output was overwritten.
