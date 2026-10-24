# ☀️ HEL1OS — Solar Flare Nowcasting (Proposal)

![Status](https://img.shields.io/badge/Status-Proposal_%2B_Research_in_progress-orange?style=for-the-badge)
![Python](https://img.shields.io/badge/Python-planned-3776AB?style=for-the-badge&logo=python&logoColor=white)
![ISRO](https://img.shields.io/badge/ISRO-Bharatiya_Antariksh_Hackathon-0B3D91?style=for-the-badge)

Our team's proposal for **ISRO's Bharatiya Antariksh Hackathon**: an ML system that **nowcasts** (detects in
real time) and **forecasts** (warns minutes ahead) solar flares from the X-ray instruments on ISRO's
**Aditya-L1** mission, with alerts for communication risk.

> **This is a proposal, not a built system.** We submitted the idea and a technical design to the
> hackathon. I'm now developing the approach into a research paper.

I led a team of three and chose this problem from 15 ISRO problem statements.

---

## 🌞 The Problem

Solar flares are sudden bursts of radiation from the Sun. Strong ones (M- and X-class) can disrupt HF radio,
GPS and satellites on Earth. Aditya-L1 sits at the Sun–Earth L1 point with two X-ray spectrometers:

| Instrument | Band | Traces |
|------------|------|--------|
| **SoLEXS** | 2–22 keV (soft X-ray) | hot flare plasma |
| **HEL1OS** | 10–150 keV (hard X-ray) | energetic electrons — often spikes *first* |

The hard-X-ray spike usually lines up with the rise of the soft-X-ray curve (the **Neupert effect**), so
combining both instruments gives an early, reliable signal that a flare is starting.

---

## 🧠 Proposed Approach

```
SoLEXS + HEL1OS light curves ─┐
NOAA active-region lists ─────┼──► features ──► Random Forest + XGBoost ensemble
SDO/HMI magnetograms ─────────┘                         │
                                    multi-instrument confidence fusion
                                                        │
                       four-tier alert, from Watch up to Emergency
```

- **Data:** Aditya-L1 Level-1 light curves (ISSDC PRADAN portal), GOES X-ray flux for flare-class labels,
  NOAA active regions and SDO/HMI magnetograms for context.
- **Model:** a Random Forest + XGBoost ensemble, with confidence fusion across the instruments.
- **Output:** a four-tier alert level for communication-risk warnings.

---

## 📁 What's in This Repo

| File | Contents |
|------|----------|
| `Aditya.docx` | Background and planning notes: flare physics, Aditya-L1 data access, tools, and the phased plan we followed while preparing the submission |

The full 31-page technical design document was submitted to the judges and isn't published here.

---

## 👩‍💻 Author

**Ashritha Abbabathula** — team lead
[![GitHub](https://img.shields.io/badge/GitHub-Ashritha--Abbabathula-181717?style=flat&logo=github)](https://github.com/Ashritha-Abbabathula)
