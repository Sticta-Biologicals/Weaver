import './App.css';
import React from 'react'
//import Massive from './Massive';
import PlasmidsTable from './PlasmidsTable';
// comment on production mode
// import testData from './test.json'

const test_mode = false
// comment on test mode
const testData = []

class App extends React.Component {
    constructor(props) {
        super(props)
        this.state = {
            data: false
        }
    }

    componentDidMount() {
        const axios = require('axios');
        const url = '/inventory/api/plasmids/'

        if (test_mode)
            this.setState({
                data: testData
            })
        else
            axios.get(url)
                .then((response) => {
                    if (response.data) {
                        if (response.data.plasmids) {
                            response.data.plasmids.sort((a, b) => a.ix - b.ix)
                        }
                        this.setState({
                            data: response.data
                        })
                    }
                })
    }

    componentDidUpdate() {
        window.onReady()
        if (!test_mode)
            window.do_filter_default()
    }

    render() {
        if (this.state.data) {
            if (this.state.data.plasmids)
                return <div id='plasmid-wrapper'>
                    <div id='plasmid-massive-wrapper'>
                        {/*<Massive plasmids={this.state.data.plasmids} />*/}
                    </div>
                    <div id='plasmid-table-wrapper'>
                        <PlasmidsTable data={this.state.data} />
                    </div>
                </div>
            else
                return <div className="alert alert-info">
                    <i className="bi bi-emoji-frown"></i> No plasmids
                </div>
        } else {
            return <div className="alert alert-info">
                <div className="spinner-grow spinner-grow-sm" role="status">
                    <span className="visually-hidden">...</span>
                </div> Loading plasmids
            </div>
        }
    }
}

export default App;
