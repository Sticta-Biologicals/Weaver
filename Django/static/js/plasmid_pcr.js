(function () {
    const selector = document.getElementById('pcr-primer-selector')
    if (!selector) return

    const form = document.getElementById('pcr-primer-selector-form')
    const filterUrl = selector.dataset.filterUrl
    const status = document.getElementById('pcr-filter-status')
    const activeFilters = document.getElementById('pcr-active-filters')
    const message = document.getElementById('pcr-primer-selector-message')
    const submit = document.getElementById('pcr-primer-submit')
    const hidden = {
        forward: document.getElementById('pcr-selected-forward'),
        reverse: document.getElementById('pcr-selected-reverse')
    }
    const customInputs = {
        forward: document.getElementById('pcr-custom-forward'),
        reverse: document.getElementById('pcr-custom-reverse')
    }
    const state = { forward: '', reverse: '' }
    const allowedIds = { forward: null, reverse: null }
    const searchTerms = { forward: '', reverse: '' }
    const options = {
        forward: Array.from(document.querySelectorAll('.pcr-primer-option[data-direction="forward"]')),
        reverse: Array.from(document.querySelectorAll('.pcr-primer-option[data-direction="reverse"]'))
    }
    let requestNumber = 0
    let filtering = false
    const filterTimers = { forward: null, reverse: null }

    function otherDirection(direction) {
        return direction === 'forward' ? 'reverse' : 'forward'
    }

    function hasCustom(direction) {
        return Boolean(customInputs[direction].value.trim())
    }

    function hasSelection(direction) {
        return Boolean(state[direction] || hasCustom(direction))
    }

    function visibleOptions(direction) {
        return options[direction].filter((option) => !option.hidden)
    }

    function applyVisibility(direction) {
        options[direction].forEach((option) => {
            const matchesPair = !allowedIds[direction] || allowedIds[direction].has(option.dataset.primerId)
            const matchesSearch = !searchTerms[direction] || option.dataset.search.toLowerCase().includes(searchTerms[direction])
            const visible = matchesPair && matchesSearch
            option.hidden = !visible
            option.classList.toggle('d-none', !visible)
            option.setAttribute('aria-disabled', visible ? 'false' : 'true')
        })
    }

    function setAllowed(direction, ids) {
        allowedIds[direction] = ids
        applyVisibility(direction)
    }

    function setActive(direction) {
        options[direction].forEach((option) => {
            const selected = option.dataset.primerId === state[direction]
            option.classList.toggle('active', selected)
            option.setAttribute('aria-checked', selected ? 'true' : 'false')
        })
    }

    function updateSubmit() {
        submit.disabled = filtering || !hasSelection('forward') || !hasSelection('reverse')
        if (hasSelection('forward') && hasSelection('reverse')) {
            message.textContent = 'Pair selected. You can design the PCR.'
        } else if (hasSelection('forward') || hasSelection('reverse')) {
            message.textContent = 'Now select a primer from the opposite panel.'
        } else {
            message.textContent = 'Select both primers to continue.'
        }
    }

    function clearSelection(direction) {
        state[direction] = ''
        hidden[direction].value = ''
        customInputs[direction].value = ''
        setActive(direction)
    }

    function selectionParameters(direction) {
        const params = 'min_size=' + encodeURIComponent(document.getElementById('pcr-min-size').value.trim()) +
            '&max_size=' + encodeURIComponent(document.getElementById('pcr-max-size').value.trim()) +
            '&max_tm_diff=' + encodeURIComponent(document.getElementById('pcr-max-tm-diff').value.trim())
        if (state[direction]) {
            return 'primer_id=' + encodeURIComponent(state[direction]) + '&' + params
        }
        return 'primer_sequence=' + encodeURIComponent(customInputs[direction].value.trim()) +
            '&direction=' + encodeURIComponent(direction) + '&' + params
    }

    function updateActiveFilters() {
        if (!activeFilters) return
        const minSize = document.getElementById('pcr-min-size').value.trim()
        const maxSize = document.getElementById('pcr-max-size').value.trim()
        const maxTmDiff = document.getElementById('pcr-max-tm-diff').value.trim()
        const productSize = minSize && maxSize ? minSize + '–' + maxSize :
            (minSize ? minSize + '+' : (maxSize ? '≤' + maxSize : 'Any'))
        const tmDiff = maxTmDiff ? maxTmDiff + ' °C' : 'Any'
        activeFilters.textContent = 'Active filters: Product size ' + productSize + ' bp · Maximum ΔTm ' + tmDiff
    }

    function updateStatus() {
        status.textContent = 'F: ' + visibleOptions('forward').length +
            ' · R: ' + visibleOptions('reverse').length
    }

    async function filterFrom(direction) {
        const targetDirection = otherDirection(direction)
        const currentRequest = ++requestNumber

        if (!hasSelection(direction)) {
            filtering = false
            if (hasSelection(targetDirection)) {
                return filterFrom(targetDirection)
            }
            setAllowed(targetDirection, null)
            updateStatus()
            updateSubmit()
            return
        }

        filtering = true
        updateSubmit()
        status.textContent = 'Searching for compatible pairs…'

        try {
            const response = await fetch(filterUrl + '?' + selectionParameters(direction), {
                headers: { 'X-Requested-With': 'XMLHttpRequest' }
            })
            const data = await response.json()
            if (currentRequest !== requestNumber) return
            if (!response.ok) throw new Error(data.error || 'Could not filter primers.')

            const allowedIds = new Set(data.compatible_ids || [])
            setAllowed(targetDirection, allowedIds)
            if (state[targetDirection] && !allowedIds.has(state[targetDirection])) {
                state[targetDirection] = ''
                hidden[targetDirection].value = ''
                setActive(targetDirection)
            }
            filtering = false
            status.textContent = allowedIds.size + ' compatible primer' + (allowedIds.size === 1 ? '' : 's')
            updateSubmit()
        } catch (error) {
            if (currentRequest !== requestNumber) return
            filtering = false
            setAllowed(targetDirection, new Set())
            status.textContent = 'Could not update the filter.'
            message.textContent = error.message
            updateSubmit()
        }
    }

    function scheduleFilter(direction) {
        clearTimeout(filterTimers[direction])
        filterTimers[direction] = setTimeout(() => filterFrom(direction), 200)
    }

    Object.keys(options).forEach((direction) => {
        options[direction].forEach((option) => {
            option.addEventListener('click', () => {
                if (state[direction] === option.dataset.primerId) {
                    clearSelection(direction)
                    filterFrom(direction)
                    return
                }
                state[direction] = option.dataset.primerId
                hidden[direction].value = state[direction]
                customInputs[direction].value = ''
                setActive(direction)
                filterFrom(direction)
            })
        })
    })

    document.querySelectorAll('.pcr-primer-search').forEach((input) => {
        input.addEventListener('input', () => {
            const direction = input.dataset.direction
            searchTerms[direction] = input.value.trim().toLowerCase()
            applyVisibility(direction)
            updateStatus()
        })
    })

    Object.keys(customInputs).forEach((direction) => {
        customInputs[direction].addEventListener('input', () => {
            if (customInputs[direction].value.trim()) {
                state[direction] = ''
                hidden[direction].value = ''
                setActive(direction)
            }
            scheduleFilter(direction)
            updateSubmit()
        })
    })

    document.querySelectorAll('.pcr-clear-selection').forEach((button) => {
        button.addEventListener('click', () => {
            const direction = button.dataset.direction
            clearSelection(direction)
            filterFrom(direction)
        })
    })

    const parametersToggle = document.getElementById('pcr-parameters-toggle')
    const parametersBody = document.getElementById('pcr-parameters-body')
    if (parametersToggle && parametersBody) {
        parametersToggle.addEventListener('click', () => {
            const collapsed = parametersBody.classList.toggle('is-collapsed')
            parametersToggle.setAttribute('aria-expanded', collapsed ? 'false' : 'true')
            const icon = parametersToggle.querySelector('i')
            if (icon) icon.className = collapsed ? 'bi bi-chevron-down' : 'bi bi-chevron-up'
        })
    }

    const applyFilters = document.getElementById('pcr-apply-filters')
    if (applyFilters) {
        applyFilters.addEventListener('click', () => {
            if (hasSelection('forward')) {
                updateActiveFilters()
                filterFrom('forward')
            } else if (hasSelection('reverse')) {
                updateActiveFilters()
                filterFrom('reverse')
            } else {
                message.textContent = 'Select a primer before applying pair filters.'
            }
        })
    }

    form.addEventListener('submit', (event) => {
        hidden.forward.value = state.forward
        hidden.reverse.value = state.reverse
        const customForward = document.querySelector('input[name="primer_f_seq"]')
        const customReverse = document.querySelector('input[name="primer_r_seq"]')
        customForward.value = customInputs.forward.value.trim()
        customReverse.value = customInputs.reverse.value.trim()
        if (!hasSelection('forward') || !hasSelection('reverse') || filtering) {
            event.preventDefault()
        }
    })

    updateStatus()
    updateSubmit()
})()
