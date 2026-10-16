#!/usr/bin/env python3
"""Graphical interface for zotu_taxonomy.py (Tkinter, no extra dependencies)."""

import os
import queue
import subprocess
import sys
import threading
import time
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import zotu_taxonomy as zt

APP_TITLE = "zOTU Taxonomy (NCBI 16S type strains)"


def open_path(path):
    if sys.platform.startswith("win"):
        os.startfile(path)  # noqa: S606 - opens the result in Excel
    elif sys.platform == "darwin":
        subprocess.Popen(["open", path])
    else:
        subprocess.Popen(["xdg-open", path])


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.minsize(680, 520)
        self.events = queue.Queue()
        self.worker = None
        self.result_path = None
        self.started = None
        self.searching = False

        self.fasta = tk.StringVar()
        self.output = tk.StringVar()
        self.top = tk.IntVar(value=zt.DEFAULTS["top"])
        self.species_id = tk.DoubleVar(value=zt.SPECIES_ID)
        self.genus_id = tk.DoubleVar(value=zt.GENUS_ID)
        self.tie_delta = tk.DoubleVar(value=0.0)
        self.min_cov = tk.DoubleVar(value=zt.DEFAULTS["min_cov"])
        self.threads = tk.IntVar(value=zt.DEFAULTS["threads"])
        self.trim = tk.BooleanVar(value=True)

        self._build()
        self._refresh_db_label()
        self.after(100, self._poll)

    # ------------------------------------------------------------------ UI
    def _build(self):
        pad = {"padx": 8, "pady": 4}
        files = ttk.LabelFrame(self, text="Αρχεία")
        files.pack(fill="x", **pad)
        ttk.Label(files, text="zOTU FASTA:").grid(row=0, column=0, sticky="w", **pad)
        ttk.Entry(files, textvariable=self.fasta).grid(row=0, column=1, sticky="ew", **pad)
        ttk.Button(files, text="Επιλογή…", command=self._pick_fasta).grid(row=0, column=2, **pad)
        ttk.Label(files, text="Αποτέλεσμα (.xlsx):").grid(row=1, column=0, sticky="w", **pad)
        ttk.Entry(files, textvariable=self.output).grid(row=1, column=1, sticky="ew", **pad)
        ttk.Button(files, text="Αποθήκευση ως…", command=self._pick_output).grid(row=1, column=2, **pad)
        files.columnconfigure(1, weight=1)

        opts = ttk.LabelFrame(self, text="Ρυθμίσεις")
        opts.pack(fill="x", **pad)
        rows = [
            ("Είδη ανά zOTU:", self.top, 1, 20, 1),
            ("Όριο είδους (%):", self.species_id, 90, 100, 0.05),
            ("Όριο γένους (%):", self.genus_id, 80, 100, 0.1),
            ("Ανοχή ισοπαλίας (%):", self.tie_delta, 0, 3, 0.1),
            ("Threads (πυρήνες):", self.threads, 1, max(os.cpu_count() or 1, 1), 1),
            ("Ελάχ. κάλυψη zOTU (0-1):", self.min_cov, 0.5, 1.0, 0.01),
        ]
        for i, (label, var, lo, hi, step) in enumerate(rows):
            r, c = divmod(i, 2)
            ttk.Label(opts, text=label).grid(row=r, column=c * 2, sticky="w", **pad)
            ttk.Spinbox(opts, textvariable=var, from_=lo, to=hi, increment=step,
                        width=8).grid(row=r, column=c * 2 + 1, sticky="w", **pad)
        ttk.Checkbutton(opts, text="Αυτόματο κόψιμο υπολειμμάτων primer",
                        variable=self.trim).grid(row=3, column=0, columnspan=2, sticky="w", **pad)

        db = ttk.Frame(self)
        db.pack(fill="x", **pad)
        self.db_label = ttk.Label(db, text="")
        self.db_label.pack(side="left")
        self.update_btn = ttk.Button(db, text="Ενημέρωση βάσης NCBI", command=self._update_db)
        self.update_btn.pack(side="right")

        run = ttk.Frame(self)
        run.pack(fill="x", **pad)
        self.run_btn = ttk.Button(run, text="▶  Εκτέλεση", command=self._run)
        self.run_btn.pack(side="left")
        self.open_btn = ttk.Button(run, text="Άνοιγμα Excel", state="disabled",
                                   command=lambda: open_path(self.result_path))
        self.open_btn.pack(side="left", padx=8)
        self.status = ttk.Label(run, text="Έτοιμο")
        self.status.pack(side="left", padx=8)

        self.bar = ttk.Progressbar(self, maximum=100)
        self.bar.pack(fill="x", **pad)

        logf = ttk.LabelFrame(self, text="Καταγραφή")
        logf.pack(fill="both", expand=True, **pad)
        self.log_box = tk.Text(logf, height=12, wrap="word", state="disabled")
        sb = ttk.Scrollbar(logf, command=self.log_box.yview)
        self.log_box.configure(yscrollcommand=sb.set)
        self.log_box.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")

    def _refresh_db_label(self):
        d = zt.db_date()
        self.db_label.config(text=f"Βάση NCBI: κατέβηκε {d}" if d else
                             "Βάση NCBI: θα κατέβει στην πρώτη εκτέλεση (~8 MB)")

    def _pick_fasta(self):
        path = filedialog.askopenfilename(
            title="Επιλογή zOTU FASTA",
            filetypes=[("FASTA", "*.fa *.fasta *.fna *.fas *.txt"), ("Όλα", "*.*")])
        if path:
            self.fasta.set(path)
            self.output.set(os.path.splitext(path)[0] + "_taxonomy.xlsx")

    def _pick_output(self):
        path = filedialog.asksaveasfilename(
            title="Αποθήκευση αποτελέσματος", defaultextension=".xlsx",
            filetypes=[("Excel", "*.xlsx")], initialfile=os.path.basename(self.output.get()))
        if path:
            self.output.set(path)

    def _log(self, text):
        self.log_box.configure(state="normal")
        self.log_box.insert("end", text + "\n")
        self.log_box.see("end")
        self.log_box.configure(state="disabled")

    def _busy(self, busy):
        state = "disabled" if busy else "normal"
        self.run_btn.config(state=state)
        self.update_btn.config(state=state)

    # ------------------------------------------------------------ actions
    def _start(self, target, *args):
        self._busy(True)
        self.open_btn.config(state="disabled")
        self.bar.config(value=0)
        self.started = time.time()
        self.worker = threading.Thread(target=target, args=args, daemon=True)
        self.worker.start()

    def _run(self):
        if self.worker and self.worker.is_alive():
            return
        fasta = self.fasta.get().strip()
        if not fasta or not os.path.isfile(fasta):
            messagebox.showerror(APP_TITLE, "Διάλεξε πρώτα ένα αρχείο FASTA.")
            return
        try:
            opts = dict(top=int(self.top.get()), species_id=float(self.species_id.get()),
                        genus_id=float(self.genus_id.get()), tie_delta=float(self.tie_delta.get()),
                        threads=int(self.threads.get()), trim="auto" if self.trim.get() else "none",
                        min_cov=float(self.min_cov.get()))
            if not 0 < opts["min_cov"] <= 1:
                raise ValueError
        except (tk.TclError, ValueError):
            messagebox.showerror(APP_TITLE, "Μη έγκυρη τιμή στις ρυθμίσεις.")
            return
        output = self.output.get().strip() or None
        self._log(f"--- Εκτέλεση: {os.path.basename(fasta)} ---")
        self._start(self._work_run, fasta, output, opts)

    def _update_db(self):
        if self.worker and self.worker.is_alive():
            return
        self._log("--- Ενημέρωση βάσης NCBI ---")
        self._start(self._work_update)

    # ------------------------------------------------- worker thread side
    def _emit(self, kind, value=None):
        self.events.put((kind, value))

    def _work_run(self, fasta, output, opts):
        try:
            path, counts = zt.run_pipeline(
                fasta, output, log=lambda t: self._emit("log", t),
                progress=lambda p: self._emit("progress", p), **opts)
            self._emit("done", (path, counts))
        except Exception as e:  # report every failure in the window
            self._emit("error", f"{type(e).__name__}: {e}")

    def _work_update(self):
        try:
            zt.ensure_db(update=True, log=lambda t: self._emit("log", t))
            self._emit("db_done")
        except Exception as e:
            self._emit("error", f"{type(e).__name__}: {e}")

    # ----------------------------------------------------- UI thread side
    def _poll(self):
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == "log":
                    self._log(value)
                    if value.startswith("Running:"):
                        # vsearch reports no percentage when not attached to a console
                        self.searching = True
                        self.bar.config(mode="indeterminate")
                        self.bar.start(15)
                    elif "downloading" in value:
                        self.status.config(text="Λήψη…")
                elif kind == "progress":
                    self._stop_bar()
                    self.bar.config(value=value)
                elif kind == "done":
                    path, counts = value
                    self.result_path = path
                    self._stop_bar(100)
                    self.status.config(text=f"Ολοκληρώθηκε σε {self._elapsed()}")
                    self.open_btn.config(state="normal")
                    self._busy(False)
                    self._refresh_db_label()
                    summary = "\n".join(f"{k}: {v}" for k, v in counts.items())
                    messagebox.showinfo(APP_TITLE, f"Ολοκληρώθηκε.\n\n{summary}\n\n{path}")
                elif kind == "db_done":
                    self.status.config(text="Η βάση ενημερώθηκε")
                    self._busy(False)
                    self._refresh_db_label()
                elif kind == "error":
                    self._stop_bar(0)
                    self._log("ΣΦΑΛΜΑ: " + value)
                    self.status.config(text="Σφάλμα")
                    self._busy(False)
                    messagebox.showerror(APP_TITLE, value)
        except queue.Empty:
            pass
        if self.searching:
            self.status.config(text=f"Αναζήτηση στη βάση…  {self._elapsed()}")
        self.after(150, self._poll)

    def _stop_bar(self, value=None):
        self.searching = False
        self.bar.stop()
        self.bar.config(mode="determinate")
        if value is not None:
            self.bar.config(value=value)

    def _elapsed(self):
        s = int(time.time() - (self.started or time.time()))
        return f"{s // 60}:{s % 60:02d}"


def selftest(fasta, output):
    """Headless run used by the build to check the packaged .exe (no window)."""
    with open(output + ".log", "w", encoding="utf-8") as log:
        try:
            zt.run_pipeline(fasta, output, log=lambda t: print(t, file=log, flush=True))
        except Exception as e:
            print(f"FAILED: {type(e).__name__}: {e}", file=log)
            sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "--selftest":
        selftest(sys.argv[2], sys.argv[3])
    else:
        App().mainloop()
