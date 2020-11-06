# Before you can use the jobs API, you need to set up an access token.
# Log in to the Quantum Experience. Under "Account", generate a personal
# access token. Replace "None" below with the quoted token string.
# Uncomment the APItoken variable, and you will be ready to go.

# this is a specifc user generated API token ... # to track stats make your own
APItoken = "f520bf68a7c88090061a75cfc7f0bcf4c2536e7d1e0662337afd34d3e30bae792008499ec8f60d61bb931dce765ad9b6619f3648802a7cfb6477745ded54261c"

config = {
  "url": 'https://quantumexperience.ng.bluemix.net/api'
}

if 'APItoken' not in locals():
  raise Exception("Please set up your access token. See Qconfig.py.")
