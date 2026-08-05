import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// https://vite.dev/config/
export default defineConfig({
  plugins: [react()],
  base: '/static/experiments-reactflow/',
  build: {
    manifest: false,
    rollupOptions: {
      output: {
        entryFileNames: 'experiments-flow.js',
        chunkFileNames: 'experiments-flow-[name].js',
        assetFileNames: 'experiments-flow.[ext]',
      },
    },
  },
})
