import './App.css';
import React from 'react'

class App extends React.Component {
    constructor(props) {
        super(props)
        this.state = {
            glycerolstocks_render: '',
            loading: true
        }
    }

    loadData() {
        const axios = require('axios');
        // const url = 'http://192.168.1.64:8000/inventory/api/glycerolstocks'
        const url = '/inventory/api/glycerolstocks'
        axios.get(url)
            .then((response) => {
                let output = []
                if (response.data.glycerolstocks) {
                    response.data.glycerolstocks.forEach((glycerolstock) => {
                        const table_filters = response.data.table_filters
                        let glycerolstock_level = ""
                        if (glycerolstock.pl != null) glycerolstock_level = " filter-l" + glycerolstock.pl
                        let glycerolstock_type = ""
                        if (glycerolstock.pt != null) glycerolstock_type = " filter-t" + glycerolstock.pt
                        let table_filters_output = ""
                        table_filters.forEach((table_filter) => {
                            if (table_filter[0] === 'startswith') {
                                table_filter[2].forEach((table_filter_op) => {
                                    if (glycerolstock.pn.toLowerCase().startsWith(table_filter_op[1])) {
                                        table_filters_output += " filter-" + table_filter_op[1]
                                    }
                                })
                            }
                        })
                        let plasmid_output = "No plasmid"
                        if (glycerolstock.pi)
                            plasmid_output =
                                <a href={"/inventory/plasmid/" + glycerolstock.pi} className="btn btn-outline-secondary"
                                   role="button"><span>{glycerolstock.pn}</span><span class="plasmid_list-id badge text-bg-light border ms-1">{glycerolstock.pix}</span></a>

                        let glycerolstock_edit_output = ""
                        if(glycerolstock.p)
                            glycerolstock_edit_output = <a href={"/inventory/glycerolstock/edit/" + glycerolstock.i}
                                                           className="btn text-secondary" role="button" data-bs-toggle="tooltip" data-bs-placement="top" title="Edit"><i
                                className="bi bi-pencil-fill"></i></a>

                        let gn_name = <span>{glycerolstock.s}</span>

                        if(glycerolstock.pn)
                            gn_name = <span>{glycerolstock.s + " / " + glycerolstock.pn}</span>
                        let plasmid_icon = ""
                        if (glycerolstock.pcs === 'v') {
                            // verified
                            plasmid_icon = <i className="bi bi-check-circle ms-2" data-bs-toggle="tooltip" data-bs-placement="top" title="Validated"></i>
                        } else if(glycerolstock.pcs === 'r') {
                            // reference
                            plasmid_icon = <i className="bi bi-bookmarks ms-2" data-bs-toggle="tooltip" data-bs-placement="top" title="Reference"></i>
                        } else if(glycerolstock.pcs === 'c') {
                            // under construction
                            plasmid_icon = <i className="bi bi-hammer ms-2" data-bs-toggle="tooltip" data-bs-placement="top" title="Under construction"></i>
                        }

                        output.push(<tr
                            class={"filter-item" + glycerolstock_level + glycerolstock_type + table_filters_output}>
                            <td>
                                <a class="btn btn-success table-search-search_on me-1" data-search-all={glycerolstock.s + glycerolstock.pn + glycerolstock.pix} data-search-name={glycerolstock.s + glycerolstock.pn} data-search-idx={glycerolstock.pix} role="button" href={"/inventory/glycerolstock/" + glycerolstock.i}>
                                    {gn_name}
                                    <span class="plasmid_list-id badge text-bg-light text-success ms-1">{glycerolstock.pix}</span>
                                </a>
                            </td>
                            <td>
                                {glycerolstock_edit_output}
                                <a href={"/inventory/glycerolstock/label/" + glycerolstock.i} class="btn text-info me-1"
                                   role="button" data-bs-toggle="tooltip" data-bs-placement="top" title="Print label"><i class="bi bi-tag-fill"></i></a>
                            </td>
                            <td>
                                {glycerolstock.s}
                            </td>
                            <td>
                                {plasmid_output}
                                {plasmid_icon}
                            </td>
                            <td>
                                {glycerolstock.br}{glycerolstock.bc}
                            </td>
                            <td>
                                {glycerolstock.bn}
                            </td>
                            <td>
                                {glycerolstock.bl}
                            </td>
                        </tr>)
                    })
                    output = <table id="glycerolstocks-table" class="table table-striped table-hover sortable table-search-target">
                        <thead>
                        <tr>
                            <th scope="col">Stock</th>
                            <th scope="col">Actions</th>
                            <th scope="col">Strain</th>
                            <th scope="col">Plasmid</th>
                            <th scope="col">Position</th>
                            <th scope="col">Box</th>
                            <th scope="col">Location</th>
                        </tr>
                        </thead>
                        <tbody>{output}</tbody>
                    </table>
                } else {
                    output = <div className="alert alert-info">
                        <i className="bi bi-emoji-frown"></i> No glycerolstocks
                    </div>
                }
                this.setState({
                    glycerolstocks_render: output,
                    loading: false
                })
                window.onReady()
            })
    }

    componentDidMount() {
        this.loadData()
    }

    render() {
        let render = <div className="alert alert-info">
            <div className="spinner-grow spinner-grow-sm" role="status">
                <span className="visually-hidden">...</span>
            </div> Loading glycerol stocks
        </div>
        if (!this.state.loading) {
            render = this.state.glycerolstocks_render
        }
        return (<div>{render}</div>)
    }
}

export default App;
