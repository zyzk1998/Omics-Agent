!macro customInstall
  nsExec::ExecToStack 'taskkill /F /IM OmicsAgent.exe'
  nsExec::ExecToStack 'taskkill /F /IM local_sidecar.exe'
  Sleep 1000
!macroend
