(function (root) {
    "use strict";

    function normalizeHeader(value) {
        var text = String(value || "").replace(/^\uFEFF/, "").trim().toLowerCase();
        if (text.normalize) {
            text = text.normalize("NFD").replace(/[\u0300-\u036f]/g, "");
        }
        return text.replace(/[^a-z0-9]+/g, " ").trim();
    }

    function splitLine(line) {
        if (line.indexOf("\t") !== -1) {
            return line.split("\t");
        }
        return line.trim().split(/\s{2,}/);
    }

    function findColumn(headers, matcher) {
        for (var index = 0; index < headers.length; index += 1) {
            if (matcher(headers[index])) {
                return index;
            }
        }
        return -1;
    }

    function parseSampleId(rawValue) {
        var value = String(rawValue || "").trim();
        var match = value.match(/^(\d+)\s*c\s*(\d+)$/i);
        if (!match) {
            match = value.match(/^(\d+)\s+(\d+)$/);
        }
        if (!match) {
            return null;
        }
        return { identifier: match[1], colony: match[2] };
    }

    function normalizeConcentration(rawValue) {
        var value = String(rawValue || "").trim().replace(/\s/g, "");
        if (!value) {
            return "";
        }
        value = value.replace(/[^0-9,\.\-+]/g, "");
        if (value.indexOf(",") !== -1 && value.indexOf(".") !== -1) {
            if (value.lastIndexOf(",") > value.lastIndexOf(".")) {
                value = value.replace(/\./g, "").replace(",", ".");
            } else {
                value = value.replace(/,/g, "");
            }
        } else {
            value = value.replace(",", ".");
        }
        var numericValue = Number(value);
        return Number.isFinite(numericValue) && numericValue > 0 ? String(numericValue) : "";
    }

    function padNumber(value) {
        return String(value).length === 1 ? "0" + value : String(value);
    }

    function normalizeDate(rawValue) {
        var value = String(rawValue || "").trim();
        var match = value.match(/^(\d{4})[-\/]([01]?\d)[-\/]([0-3]?\d)/);
        if (!match) {
            match = value.match(/^([0-3]?\d)[-\/]([01]?\d)[-\/](\d{4})/);
            if (match) {
                match = [match[0], match[3], match[2], match[1]];
            }
        }
        if (!match) {
            return "";
        }

        var year = Number(match[1]);
        var month = Number(match[2]);
        var day = Number(match[3]);
        var parsed = new Date(Date.UTC(year, month - 1, day));
        if (
            !Number.isFinite(year) ||
            !Number.isFinite(month) ||
            !Number.isFinite(day) ||
            parsed.getUTCFullYear() !== year ||
            parsed.getUTCMonth() !== month - 1 ||
            parsed.getUTCDate() !== day
        ) {
            return "";
        }
        return year + "-" + padNumber(month) + "-" + padNumber(day);
    }

    function parseNanodropText(rawText) {
        var text = String(rawText || "").replace(/^\uFEFF/, "").replace(/\r/g, "");
        var lines = text.split("\n").filter(function (line) {
            return line.trim();
        });
        if (!lines.length) {
            return { rows: [], invalidSamples: [], errors: ["No NanoDrop data was found."] };
        }

        var headerLineIndex = -1;
        var headers = [];
        for (var lineIndex = 0; lineIndex < lines.length; lineIndex += 1) {
            var candidateHeaders = splitLine(lines[lineIndex]).map(normalizeHeader);
            if (findColumn(candidateHeaders, function (header) {
                return header === "sample id" || header === "sampleid";
            }) !== -1) {
                headerLineIndex = lineIndex;
                headers = candidateHeaders;
                break;
            }
        }
        if (headerLineIndex === -1) {
            return {
                rows: [],
                invalidSamples: [],
                errors: ["The pasted table must include a Sample ID column."]
            };
        }

        var sampleColumn = findColumn(headers, function (header) {
            return header === "sample id" || header === "sampleid";
        });
        var concentrationColumn = findColumn(headers, function (header) {
            return header.indexOf("nucleic acid") !== -1 ||
                header.indexOf("concentration") !== -1 ||
                /^conc(?:$| )/.test(header);
        });
        var dateColumn = findColumn(headers, function (header) {
            return header === "date" || header.indexOf("date and time") !== -1 || header.indexOf("date time") !== -1;
        });
        var missingColumns = [];
        if (concentrationColumn === -1) {
            missingColumns.push("Nucleic Acid/concentration");
        }
        if (dateColumn === -1) {
            missingColumns.push("Date and Time/date");
        }
        if (missingColumns.length) {
            return {
                rows: [],
                invalidSamples: [],
                errors: ["Missing column(s): " + missingColumns.join(", ") + "."]
            };
        }

        var rows = [];
        var invalidSamples = [];
        for (var dataIndex = headerLineIndex + 1; dataIndex < lines.length; dataIndex += 1) {
            var cells = splitLine(lines[dataIndex]);
            var rawSampleId = cells[sampleColumn] || "";
            if (!rawSampleId.trim() || /^#?$/.test(rawSampleId.trim())) {
                continue;
            }
            var sample = parseSampleId(rawSampleId);
            if (!sample) {
                invalidSamples.push(rawSampleId.trim());
                continue;
            }
            rows.push({
                identifier: sample.identifier,
                colony: sample.colony,
                concentration: normalizeConcentration(cells[concentrationColumn]),
                date: normalizeDate(cells[dateColumn])
            });
        }

        return { rows: rows, invalidSamples: invalidSamples, errors: [] };
    }

    function setValue(row, selector, value) {
        var element = row.querySelector(selector);
        if (element) {
            element.value = value || "";
        }
    }

    function initBatchPrints() {
        if (!root.document) {
            return;
        }
        var page = root.document.querySelector(".batch-prints-page");
        var rowsContainer = root.document.querySelector("#batch-print-rows");
        if (!page || !rowsContainer) {
            return;
        }

        var today = page.getAttribute("data-today") || "";
        var initialRow = rowsContainer.querySelector(".batch-print-row");
        if (!initialRow) {
            return;
        }
        var rowTemplate = initialRow.cloneNode(true);
        var status = root.document.querySelector("#batch-nanodrop-status");
        var pasteButton = root.document.querySelector("#batch-paste-nanodrop");
        var pasteAttempt = 0;
        var clipboardTarget = root.document.createElement("textarea");
        clipboardTarget.className = "batch-nanodrop-clipboard-target";
        clipboardTarget.setAttribute("aria-label", "Paste NanoDrop data");
        clipboardTarget.setAttribute("tabindex", "0");
        clipboardTarget.spellcheck = false;
        page.appendChild(clipboardTarget);

        function refreshConcentrationState(row) {
            var type = row.querySelector(".batch-type");
            var concentration = row.querySelector(".batch-concentration");
            if (!type || !concentration) {
                return;
            }
            if (type.value === "glycerolstocks") {
                concentration.value = "";
                concentration.readOnly = true;
                concentration.placeholder = "N/A";
            } else {
                concentration.readOnly = false;
                concentration.placeholder = "ng/ul";
            }
        }

        function clearRow(row, labelType) {
            setValue(row, ".batch-type", labelType || "plasmids");
            setValue(row, ".batch-identifier", "");
            setValue(row, ".batch-colony", "");
            setValue(row, ".batch-date", today);
            setValue(row, ".batch-concentration", "");
            refreshConcentrationState(row);
        }

        function bindRow(row) {
            refreshConcentrationState(row);
            var type = row.querySelector(".batch-type");
            var removeButton = row.querySelector(".batch-remove-row");
            if (type) {
                type.addEventListener("change", function () {
                    refreshConcentrationState(row);
                });
            }
            if (removeButton) {
                removeButton.addEventListener("click", function () {
                    if (rowsContainer.querySelectorAll(".batch-print-row").length > 1) {
                        row.remove();
                    } else {
                        clearRow(row, "plasmids");
                    }
                });
            }
        }

        function appendRow(values) {
            var row = rowTemplate.cloneNode(true);
            clearRow(row, "plasmids");
            if (values) {
                setValue(row, ".batch-identifier", values.identifier);
                setValue(row, ".batch-colony", values.colony);
                setValue(row, ".batch-concentration", values.concentration);
                setValue(row, ".batch-date", values.date);
            }
            rowsContainer.appendChild(row);
            bindRow(row);
            return row;
        }

        function replaceRows(values) {
            rowsContainer.innerHTML = "";
            values.forEach(appendRow);
            if (!values.length) {
                appendRow();
            }
        }

        function showStatus(message, kind) {
            if (!status) {
                return;
            }
            status.textContent = message;
            status.hidden = !message;
            status.className = "batch-nanodrop-status" + (kind ? " " + kind : "");
        }

        function setPasteMode(active) {
            if (!pasteButton) {
                return;
            }
            pasteButton.classList.toggle("is-listening", active);
            pasteButton.innerHTML = active ?
                '<i class="bi bi-clipboard-check"></i> Paste now (Ctrl+V)' :
                '<i class="bi bi-clipboard-plus"></i> Paste NanoDrop data';
        }

        function importText(text) {

            var result = parseNanodropText(text);
            if (!result.rows.length) {
                showStatus(result.errors.join(" ") || "No valid NanoDrop samples were found.", "is-warning");
                return;
            }
            replaceRows(result.rows);
            setPasteMode(false);
            var message = result.rows.length + " sample row" + (result.rows.length === 1 ? "" : "s") + " loaded.";
            if (result.invalidSamples.length) {
                message += " Skipped unrecognized Sample ID(s): " + result.invalidSamples.join(", ") + ".";
            }
            showStatus(message, result.invalidSamples.length ? "is-warning" : "is-success");
            rowsContainer.querySelector(".batch-identifier").focus();
        }

        rowsContainer.querySelectorAll(".batch-print-row").forEach(bindRow);
        root.document.querySelectorAll(".batch-add-entry").forEach(function (button) {
            button.addEventListener("click", function () {
                var row = appendRow();
                setValue(row, ".batch-type", button.getAttribute("data-label-type"));
                refreshConcentrationState(row);
                row.querySelector(".batch-identifier").focus();
            });
        });

        clipboardTarget.addEventListener("paste", function (event) {
            var pastedText = event.clipboardData && event.clipboardData.getData("text");
            if (pastedText) {
                event.preventDefault();
                importText(pastedText);
            }
        });

        if (pasteButton) {
            pasteButton.addEventListener("click", function () {
                setPasteMode(true);
                clipboardTarget.value = "";
                clipboardTarget.focus();
                clipboardTarget.select();
                showStatus("Ready to receive NanoDrop data. Press Ctrl+V.", "is-ready");
                if (!root.navigator || !root.navigator.clipboard || !root.navigator.clipboard.readText) {
                    return;
                }
                root.navigator.clipboard.readText().then(function (text) {
                    if (text && text.trim()) {
                        clipboardTarget.value = text;
                        importText(text);
                    }
                }).catch(function () {
                    showStatus("Ready to receive NanoDrop data. Press Ctrl+V.", "is-ready");
                });
            });
        }
    }

    var api = { parseNanodropText: parseNanodropText, parseSampleId: parseSampleId };
    root.BatchPrints = api;
    if (typeof module !== "undefined" && module.exports) {
        module.exports = api;
    }
    if (root.document) {
        if (root.document.readyState === "loading") {
            root.document.addEventListener("DOMContentLoaded", initBatchPrints);
        } else {
            initBatchPrints();
        }
    }
}(typeof window !== "undefined" ? window : globalThis));
