#!/usr/bin/env python3
"""
write_placeholder_plot.py
Writes a simple matplotlib figure with a centred message to sex_plot.png.

Usage: python3 write_placeholder_plot.py "Your message here"
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

msg = sys.argv[1] if len(sys.argv) > 1 else "Sex check unavailable"

fig, ax = plt.subplots(figsize=(8, 4))
ax.text(0.5, 0.5, msg, ha="center", va="center",
        fontsize=13, transform=ax.transAxes, wrap=True)
ax.axis("off")
plt.tight_layout()
plt.savefig("sex_plot.png", dpi=150, bbox_inches="tight")
plt.close()
